"""
    A collection of classes to make it easier to use/read output files
    generated by Starlight (from Cid Fernandes).
"""

import os
import re
from typing import Any, Type
from dataclasses import dataclass, field
import numpy as np


@dataclass
class StarlightOutput():
    """
        Data containing starlight output results (for starlight v04).
    """

    arq_obs: str = ""
    arq_base: str = ""
    arq_masks: str = ""
    arq_config: str = ""
    n_base: int = np.nan
    n_yav_components: int = np.nan
    i_fit_power_law: int = np.nan
    alpha_power_law: float = np.nan
    red_law_option: str = ""
    q_norm: float = np.nan

    l_ini: float = np.nan
    l_fin: float = np.nan
    delta_l: float = np.nan

    l_norm: float = np.nan
    llow_norm: float = np.nan
    lupp_norm: float = np.nan
    fobs_norm: float = np.nan

    llow_sn: float = np.nan
    lupp_sn: float = np.nan
    sn_in_sn_window: float = np.nan
    sn_in_norm_window: float = np.nan
    sn_err_in_sn_window: float = np.nan
    sn_err_in_norm_window: float = np.nan
    fscale_chi2: float = np.nan

    idum_orig: int = np.nan
    n0l_eff: int = np.nan
    nl_eff: int = np.nan
    ntot_cliped: int = np.nan
    clip_method: str = ""
    n_global_steps: int = np.nan
    n_chains: int = np.nan
    n_ex0s_base: int = np.nan
    clip_bug: int = np.nan
    rc_crash: int = np.nan
    burn_in: int = np.nan
    n_censored_weights: int = np.nan
    wei_nsig_threshold: float = np.nan
    wei_limit: float = np.nan
    idt_all: int = np.nan
    wdt_tot_time: float = np.nan
    wdt_usr_time: float = np.nan
    wdt_sys_time: float = np.nan

    chi2_nl_eff_ratio: float = np.nan
    adev: float = np.nan
    sum_of_x: float = np.nan
    flux_tot: float = np.nan
    m_ini_tot: float = np.nan
    m_cor_tot: float = np.nan

    v0_min: float = np.nan
    vd_min: float = np.nan
    av_min: float = np.nan
    yav_min: float = np.nan

    x_j: np.ndarray = field(default_factory=lambda:
                            np.zeros((0), dtype=float))
    m_ini_j: np.ndarray = field(default_factory=lambda:
                                np.zeros((0), dtype=float))
    m_cor_j: np.ndarray = field(default_factory=lambda:
                                np.zeros((0), dtype=float))
    age_j: np.ndarray = field(default_factory=lambda:
                              np.zeros((0), dtype=float))
    z_j: np.ndarray = field(default_factory=lambda:
                            np.zeros((0), dtype=float))
    l_m_ratio_j: np.ndarray = field(default_factory=lambda:
                                    np.zeros((0), dtype=float))
    yav_j: np.ndarray = field(default_factory=lambda:
                              np.zeros((0), dtype=float))
    m_stars: np.ndarray = field(default_factory=lambda:
                                np.zeros((0), dtype=float))
    component_j: np.ndarray = field(default_factory=lambda:
                                    np.zeros((0), dtype=float))
    a_fe_ratio_j: np.ndarray = field(default_factory=lambda:
                                     np.zeros((0), dtype=float))
    ssp_chi2r_j: np.ndarray = field(default_factory=lambda:
                                    np.zeros((0), dtype=float))
    ssp_adev_j: np.ndarray = field(default_factory=lambda:
                                   np.zeros((0), dtype=float))
    ssp_av_j: np.ndarray = field(default_factory=lambda:
                                 np.zeros((0), dtype=float))
    ssp_x_j: np.ndarray = field(default_factory=lambda:
                                np.zeros((0), dtype=float))

    l_obs: np.ndarray = field(default_factory=lambda:
                              np.zeros((0), dtype=float))
    f_obs: np.ndarray = field(default_factory=lambda:
                              np.zeros((0), dtype=float))
    f_syn: np.ndarray = field(default_factory=lambda:
                              np.zeros((0), dtype=float))
    wei: np.ndarray = field(default_factory=lambda:
                            np.zeros((0), dtype=float))


class StarlightOutputReader():
    """
        Get starlight output parameters from file (for starlight v04).
    """

    def __init__(self):
        self._re_brackets = re.compile(r"\[.*?\]")
        self._re_spaces = re.compile(r"\s+")
        self._re_trim = re.compile(r"^[\s]*|[\s]*$")
        self._file_lines = list()
        self._current_file = None

    def get_parameters(self, sl_output: str) -> StarlightOutput | None:
        """
            Read starlight output file.
        """

        out_par = StarlightOutput()

        if not self._is_file_valid(sl_output):
            return out_par

        with open(sl_output, "r", encoding="utf8") as sl_file:
            self._file_lines = sl_file.readlines()

        self._current_file = sl_output

        self._feed_simple_output(out_par)
        self._feed_array_output(out_par)

        self._clear_temp_data()

        return out_par

    def _clear_temp_data(self):
        self._file_lines.clear()
        self._current_file = None

    def _is_file_valid(self, sl_output: str) -> bool:
        return os.path.exists(sl_output) and not os.path.isdir(sl_output)

    def _feed_simple_output(self, out_par: StarlightOutput) -> None:
        out_par.arq_obs = self._get_value("[arq_obs]", str)
        out_par.arq_base = self._get_value("[arq_base]", str)
        out_par.arq_masks = self._get_value("[arq_masks]", str)
        out_par.arq_config = self._get_value("[arq_config]", str)
        out_par.n_base = self._get_value("[N_base]", int)
        out_par.n_yav_components = self._get_value(
            "[N_YAV_components = # of components with extra extinction!]", int)
        out_par.i_fit_power_law = self._get_value(
            "[i_FitPowerLaw (1/0 = Yes/No)]", int)
        out_par.alpha_power_law = self._get_value("[alpha_PowerLaw]", float)
        out_par.red_law_option = self._get_value("[red_law_option]", str)
        out_par.q_norm = self._get_value("[q_norm = A(l_norm)/A(V)]", float)

        out_par.l_ini = self._get_value("[l_ini (A)]", float)
        out_par.l_fin = self._get_value("[l_fin (A)]", float)
        out_par.delta_l = self._get_value("[dl    (A)]", float)

        out_par.l_norm = self._get_value("[l_norm (A) - for base]", float)
        out_par.llow_norm = self._get_value(
            "[llow_norm (A) - window for f_obs]", float)
        out_par.lupp_norm = self._get_value(
            "[lupp_norm (A) - window for f_obs]", float)
        out_par.fobs_norm = self._get_value(
            "[fobs_norm (in input units)]", float)

        out_par.llow_sn = self._get_value(
            "[llow_SN (A) - window for S/N]", float)
        out_par.lupp_sn = self._get_value(
            "[lupp_SN (A) - window for S/N]", float)
        out_par.sn_in_sn_window = self._get_value(
            "[S/N in S/N window]", float)
        out_par.sn_in_norm_window = self._get_value(
            "[S/N in norm. window]", float)
        out_par.sn_err_in_sn_window = self._get_value(
            "[S/N_err in S/N window]", float)
        out_par.sn_err_in_norm_window = self._get_value(
            "[S/N_err in norm. window]", float)
        out_par.fscale_chi2 = self._get_value("[fscale_chi2]", float)

        out_par.idum_orig = self._get_value("[idum_orig]", int)
        out_par.n0l_eff = self._get_value("[NOl_eff]", int)
        out_par.nl_eff = self._get_value("[Nl_eff]", int)
        out_par.ntot_cliped = self._get_value(
            "[Ntot_cliped & clip_method]", int, 0)
        out_par.clip_method = self._get_value(
            "[Ntot_cliped & clip_method]", str, 1)
        out_par.n_global_steps = self._get_value("[Nglobal_steps]", int)
        out_par.n_chains = self._get_value("[N_chains]", int)
        out_par.n_ex0s_base = self._get_value(
            "[NEX0s_base = N_base in EX0s-fits]", int)
        out_par.clip_bug = self._get_value(
            "[Clip-Bug, RC-Crash & Burn-In warning-flags, n_censored_weights, wei_nsig_threshold & wei_limit]", int, 0)
        out_par.rc_crash = self._get_value(
            "[Clip-Bug, RC-Crash & Burn-In warning-flags, n_censored_weights, wei_nsig_threshold & wei_limit]", int, 1)
        out_par.burn_in = self._get_value(
            "[Clip-Bug, RC-Crash & Burn-In warning-flags, n_censored_weights, wei_nsig_threshold & wei_limit]", int, 2)
        out_par.n_censored_weights = self._get_value(
            "[Clip-Bug, RC-Crash & Burn-In warning-flags, n_censored_weights, wei_nsig_threshold & wei_limit]", int, 3)
        out_par.wei_nsig_threshold = self._get_value(
            "[Clip-Bug, RC-Crash & Burn-In warning-flags, n_censored_weights, wei_nsig_threshold & wei_limit]", float, 4)
        out_par.wei_limit = self._get_value(
            "[Clip-Bug, RC-Crash & Burn-In warning-flags, n_censored_weights, wei_nsig_threshold & wei_limit]", float, 5)
        out_par.idt_all = self._get_value(
            "[idt_all, wdt_TotTime, wdt_UsrTime & wdt_SysTime (sec)]", int, 0)
        out_par.wdt_tot_time = self._get_value(
            "[idt_all, wdt_TotTime, wdt_UsrTime & wdt_SysTime (sec)]", float, 1)
        out_par.wdt_usr_time = self._get_value(
            "[idt_all, wdt_TotTime, wdt_UsrTime & wdt_SysTime (sec)]", float, 2)
        out_par.wdt_sys_time = self._get_value(
            "[idt_all, wdt_TotTime, wdt_UsrTime & wdt_SysTime (sec)]", float, 3)

        out_par.chi2_nl_eff_ratio = self._get_value("[chi2/Nl_eff]", float)
        out_par.adev = self._get_value("[adev (%)]", float)
        out_par.sum_of_x = self._get_value("[sum-of-x (%)]", float)
        out_par.flux_tot = self._get_value(
            "[Flux_tot (units of input spectrum!)]", float)
        out_par.m_ini_tot = self._get_value("[Mini_tot (???)]", float)
        out_par.m_cor_tot = self._get_value("[Mcor_tot (???)]", float)

        out_par.v0_min = self._get_value("[v0_min  (km/s)]", float)
        out_par.vd_min = self._get_value("[vd_min  (km/s)]", float)
        out_par.av_min = self._get_value("[AV_min  (mag)]", float)
        out_par.yav_min = self._get_value("[YAV_min (mag)]", float)

    def _feed_array_output(self, out_par: StarlightOutput) -> None:

        search_string = "# j     x_j(%)      Mini_j(%)     Mcor_j(%)     age_j(yr)     Z_j      (L/M)_j    YAV?  Mstars   component_j        a/Fe...       SSP_chi2r SSP_adev(%)   SSP_AV   SSP_x(%)"
        dtype_list = [int] + [float]*6 + [int, float, str] + [float]*5
        column_arrays = self._get_columns(search_string, dtype_list)

        out_par.x_j = column_arrays[1]
        out_par.m_ini_j = column_arrays[2]
        out_par.m_cor_j = column_arrays[3]
        out_par.age_j = column_arrays[4]
        out_par.z_j = column_arrays[5]
        out_par.l_m_ratio_j = column_arrays[6]
        out_par.yav_j = column_arrays[7]
        out_par.m_stars = column_arrays[8]
        out_par.component_j = column_arrays[9]
        out_par.a_fe_ratio_j = column_arrays[10]
        out_par.ssp_chi2r_j = column_arrays[11]
        out_par.ssp_adev_j = column_arrays[12]
        out_par.ssp_av_j = column_arrays[13]
        out_par.ssp_x_j = column_arrays[14]

        search_string = "## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei"
        dtype_list = [float] * 4
        column_arrays = self._get_columns(search_string, dtype_list, 1)

        out_par.l_obs = column_arrays[0]
        out_par.f_obs = column_arrays[1]
        out_par.f_syn = column_arrays[2]
        out_par.wei = column_arrays[3]

    def _get_columns(self, st_search: str, dtypes: list[Type], n_skip: int = 0) -> list[np.ndarray]:

        current_index = 0
        for index, line in enumerate(self._file_lines):
            if st_search in line:
                current_index = index + n_skip + 1
                break

        line_values = []

        while self._is_line_valid(current_index):
            line = self._file_lines[current_index]
            value_str = re.sub(self._re_trim, "", line)
            values = re.split(self._re_spaces, value_str)
            line_values.append(values)
            current_index += 1

        line_values = np.array(line_values, dtype=str)

        column_arrays = []
        for index, dtype in enumerate(dtypes):
            column = np.array(line_values[:, index], dtype=dtype)
            column_arrays.append(column)

        return column_arrays

    def _is_line_valid(self, index: int) -> bool:

        if index >= len(self._file_lines):
            return False

        line = self._file_lines[index]
        filtered_line = re.sub(self._re_spaces, "", line)

        return len(filtered_line) > 0

    def _get_value(self, st_search: str, dtype: Type, index: int = 0) -> Any:

        for line in self._file_lines:
            if st_search in line:
                value_str = re.sub(self._re_brackets, "", line)
                value_str = re.sub(self._re_trim, "", value_str)
                values = re.split(self._re_spaces, value_str)
                return dtype(values[index])

        message = f"Missing output compontent for {st_search}, file {self._current_file}"
        raise ValueError(message)
