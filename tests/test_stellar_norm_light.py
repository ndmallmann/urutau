"""
    Tests for the "stellar norm light" option (light normalization of
    Mage_L and MZ_L) of the STARLIGHT wrapper.

    Run from the repository root with:
        python -m unittest discover -s tests -v
"""

import unittest

import numpy as np
from astropy.io import fits

from urutau.starlight_utils import StarlightGeneric


class FakeOutput:
    """Minimal stand-in for StarlightOutput (missing scalars default to 0)."""

    def __init__(self, fc_age: float):
        # 3 SSPs (log t = 9, Z = 0.02) holding 90% of the light, plus one
        # featureless continuum holding 10%. fc_age = 0 or 1 yr, as in the
        # different STARLIGHT bases in use.
        self.invalid_file = False
        self.component_j = np.array(["ssp_a", "ssp_b", "ssp_c", "agn_fc_0.5"])
        self.x_j = np.array([30., 30., 30., 10.])
        self.age_j = np.array([1e9, 1e9, 1e9, fc_age])
        self.z_j = np.array([0.02, 0.02, 0.02, 0.])
        self.m_cor_j = np.array([40., 30., 30., 0.])
        self.m_ini_j = np.array([40., 30., 30., 0.])

    def __getattr__(self, name):
        return 0.


def popbins_header_and_data(stellar_norm_light: bool, fc_age: float):
    wrapper = StarlightGeneric.__new__(StarlightGeneric)
    wrapper._stellar_norm_light = stellar_norm_light
    wrapper._pop_age = {"xy": [0, 1e8], "xo": [1e8, 2e10]}
    wrapper._sfr_age = {}
    wrapper._fc_par = {"fc": [0., 2.]}
    wrapper._bb_par = {}
    wrapper._ret_mass_age = {}
    wrapper._galaxy_distance = 0.
    wrapper._card_index = 0
    wrapper._new_cards = list()
    wrapper._new_data = list()
    wrapper._hdu_list = fits.HDUList()
    wrapper._y_size, wrapper._x_size = 1, 1
    wrapper._matrix = np.empty((1, 1), dtype=object)
    wrapper._matrix[0, 0] = FakeOutput(fc_age)

    wrapper._add_popbins_hdu()
    header = wrapper._hdu_list[-1].header

    values = {}
    comments = {}
    for index, data in enumerate(wrapper._new_data):
        key = f"DATA{index}"
        values[header[key]] = float(data[0, 0])
        comments[header[key]] = header.comments[key]
    return header, values, comments


class TestStellarNormLight(unittest.TestCase):

    def test_total_and_stellar_light(self):
        for fc_age in (0., 1.):
            with self.subTest(fc_age=fc_age):
                h_tot, v_tot, c_tot = popbins_header_and_data(False, fc_age)
                h_ste, v_ste, c_ste = popbins_header_and_data(True, fc_age)

                # Default (total light): FC counts as log t = 0, Z = 0.
                self.assertAlmostEqual(v_tot["Mage_L"], 8.1)
                self.assertAlmostEqual(v_tot["MZ_L"], 0.018)

                # Stellar light only.
                self.assertAlmostEqual(v_ste["Mage_L"], 9.0)
                self.assertAlmostEqual(v_ste["MZ_L"], 0.02)

                # value_stellar = value_total / (1 - FC/100)
                fc = v_tot["fc"]
                self.assertAlmostEqual(fc, 10.)
                self.assertAlmostEqual(v_ste["Mage_L"], v_tot["Mage_L"] / (1. - fc / 100.))
                self.assertAlmostEqual(v_ste["MZ_L"], v_tot["MZ_L"] / (1. - fc / 100.))

                # Everything else is identical (Mage_M, MZ_M, population
                # vectors, FC fraction, ...).
                self.assertEqual(set(v_tot), set(v_ste))
                for name in v_tot:
                    if name in ("Mage_L", "MZ_L"):
                        continue
                    np.testing.assert_equal(v_tot[name], v_ste[name], err_msg=name)
                self.assertAlmostEqual(v_tot["Mage_M"], 9.0)
                self.assertAlmostEqual(v_tot["MZ_M"], 0.02)
                self.assertAlmostEqual(v_tot["xo_light"], 90.)
                self.assertAlmostEqual(v_tot["xo_mass"], 100.)

                # Traceability in the PopBins header.
                self.assertEqual(h_tot["LNORM"], "TOTAL")
                self.assertEqual(h_ste["LNORM"], "STELLAR")
                self.assertIn("total light", c_tot["Mage_L"])
                self.assertIn("stellar light", c_ste["MZ_L"])

    def test_default_is_total_light(self):
        wrapper = StarlightGeneric.__new__(StarlightGeneric)
        wrapper._stellar_norm_light = False
        self.assertAlmostEqual(
            wrapper._mean_by_light(FakeOutput(0.), np.log10(FakeOutput(0.).age_j)), 8.1)


class TestModuleDefault(unittest.TestCase):

    def test_module_default_is_off(self):
        from urutau.modules._starlight import StarlightOnUrutau
        module = StarlightOnUrutau()
        self.assertIs(module.default_parameters["stellar norm light"], False)


if __name__ == "__main__":
    unittest.main()
