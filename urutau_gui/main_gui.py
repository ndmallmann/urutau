"""
main_gui.py — Main Urutau GUI window.

A single-page interface to configure and run an Urutau pipeline
(resolution degradation, spectral/spatial resampling, dereddening,
signal-to-noise masking and STARLIGHT fitting) without editing a
run_*.py script by hand.
"""

import os
import sys
import json
import io
import traceback

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QScrollArea, QGroupBox, QCheckBox,
    QSplitter, QProgressBar, QFileDialog, QMessageBox,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont

from .constants import STYLESHEET, ACCENT, CARD_BG, MUTED, BORDER_COLOR, SUCCESS_COLOR
from .custom_widgets import (
    FilePickerRow, LabeledRow, OptionalDoubleRow, KeyRangeTable,
    TargetsTableWidget, LogConsole, CheckableGroupBox,
    labeled_double, labeled_int, labeled_combo, labeled_text,
)
from .urutau_wrapper import (
    default_config, run_pipeline, DEREDDENING_LAWS, SN_MASK_METHODS, STAT_TYPE_MAP,
)
from .base_grid_utils import read_base_components, find_agn_bin_conflicts, BaseGridError
from .script_export import generate_script


# ---------------------------------------------------------------------------
# Worker Thread — runs the Urutau pipeline in the background
# ---------------------------------------------------------------------------

class UrutauWorker(QThread):
    log_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(bool)  # True = success

    def __init__(self, cfg):
        super().__init__()
        self.cfg = cfg

    def run(self):
        old_stdout = sys.stdout
        sys.stdout = _StreamRedirector(self.log_signal)
        try:
            run_pipeline(self.cfg)
            self.finished_signal.emit(True)
        except Exception:
            tb = traceback.format_exc()
            self.log_signal.emit(f"\n❌ ERROR:\n{tb}")
            self.finished_signal.emit(False)
        finally:
            sys.stdout = old_stdout


class _StreamRedirector(io.TextIOBase):
    def __init__(self, signal):
        super().__init__()
        self._signal = signal

    def write(self, text):
        if text and text.strip():
            self._signal.emit(text)
        return len(text) if text else 0

    def flush(self):
        pass


# ---------------------------------------------------------------------------
# Main Window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    """Urutau — Main application window."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Urutau — Pipeline Configurator")
        self.resize(1240, 900)
        self.setStyleSheet(STYLESHEET)

        self.worker = None
        self._init_ui()
        self._apply_config(default_config())

    # -----------------------------------------------------------------
    # UI Construction
    # -----------------------------------------------------------------

    def _init_ui(self):
        central = QWidget(self)
        self.setCentralWidget(central)
        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        sidebar = self._create_sidebar()
        main_layout.addWidget(sidebar)

        right_side = QWidget()
        right_layout = QVBoxLayout(right_side)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(0)

        splitter = QSplitter(Qt.Vertical)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_content = QWidget()
        self.panels_layout = QVBoxLayout(scroll_content)
        self.panels_layout.setContentsMargins(20, 16, 20, 16)
        self.panels_layout.setSpacing(14)

        self._create_panel_targets_output()
        self._create_panel_input_hdus()
        self._create_panel_spatial()
        self._create_panel_resolution()
        self._create_panel_spectral_bin()
        self._create_panel_dereddening()
        self._create_panel_sn_mask()
        self._create_panel_starlight()

        self.panels_layout.addStretch()
        scroll.setWidget(scroll_content)
        splitter.addWidget(scroll)

        bottom = self._create_bottom_panel()
        splitter.addWidget(bottom)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        right_layout.addWidget(splitter)
        main_layout.addWidget(right_side, 1)

        self.statusBar().showMessage("Ready — configure the pipeline and click Run Urutau.")

    def _create_sidebar(self):
        sidebar = QWidget()
        sidebar.setFixedWidth(230)
        sidebar.setStyleSheet(
            f"background-color: {CARD_BG}; border-right: 1px solid {BORDER_COLOR};"
        )
        layout = QVBoxLayout(sidebar)
        layout.setContentsMargins(12, 20, 12, 20)
        layout.setSpacing(6)

        title = QLabel("URUTAU")
        title.setStyleSheet(f"font-size: 20px; font-weight: 800; color: {ACCENT};")
        sub = QLabel("Pipeline Configurator (prototype)")
        sub.setStyleSheet(f"font-size: 11px; color: {MUTED}; margin-bottom: 8px;")
        sub.setWordWrap(True)
        layout.addWidget(title)
        layout.addWidget(sub)

        self.nav_buttons = []
        sections = [
            ("1. Targets && Output", 0),
            ("2. Input HDUs", 1),
            ("3. Spatial Resampling", 2),
            ("4. Resolution", 3),
            ("5. Spectral Binning", 4),
            ("6. Dereddening", 5),
            ("7. S/N Mask", 6),
            ("8. Starlight (STARLIGHT)", 7),
        ]
        for text, idx in sections:
            btn = QPushButton(text)
            btn.setObjectName("navBtn")
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(lambda _, i=idx: self._scroll_to_section(i))
            self.nav_buttons.append(btn)
            layout.addWidget(btn)

        layout.addSpacing(14)

        self.btn_run = QPushButton("▶  Run Urutau")
        self.btn_run.setStyleSheet(
            f"background-color: {SUCCESS_COLOR}; color: white; font-weight: bold;"
            f"font-size: 14px; padding: 10px 16px; border-radius: 8px;"
        )
        self.btn_run.setCursor(Qt.PointingHandCursor)
        self.btn_run.clicked.connect(self._on_run)
        layout.addWidget(self.btn_run)

        layout.addSpacing(10)

        lbl_cfg = QLabel("Configuration")
        lbl_cfg.setStyleSheet(f"color: {ACCENT}; font-size: 11px; font-weight: bold;")
        layout.addWidget(lbl_cfg)

        btn_load = QPushButton("Load Config")
        btn_load.setObjectName("minorBtn")
        btn_load.setCursor(Qt.PointingHandCursor)
        btn_load.clicked.connect(self._on_load_config)
        layout.addWidget(btn_load)

        btn_save = QPushButton("Save Config")
        btn_save.setObjectName("minorBtn")
        btn_save.setCursor(Qt.PointingHandCursor)
        btn_save.clicked.connect(self._on_save_config)
        layout.addWidget(btn_save)

        btn_export_script = QPushButton("Export Script (.py)")
        btn_export_script.setObjectName("minorBtn")
        btn_export_script.setCursor(Qt.PointingHandCursor)
        btn_export_script.clicked.connect(self._on_export_script)
        layout.addWidget(btn_export_script)

        btn_reset = QPushButton("Reset to Defaults")
        btn_reset.setObjectName("minorBtn")
        btn_reset.setCursor(Qt.PointingHandCursor)
        btn_reset.clicked.connect(lambda: self._apply_config(default_config()))
        layout.addWidget(btn_reset)

        layout.addStretch()

        footer = QLabel("SAFFARY group\nUFRGS / Depto Astronomia")
        footer.setStyleSheet(f"color: {MUTED}; font-size: 11px; line-height: 1.4;")
        layout.addWidget(footer)

        return sidebar

    def _scroll_to_section(self, idx):
        targets = [
            self.grp_targets, self.grp_input, self.grp_spatial, self.grp_resolution,
            self.grp_spectral_bin, self.grp_dereddening, self.grp_sn_mask, self.grp_starlight,
        ]
        if 0 <= idx < len(targets):
            scroll = self.centralWidget().findChild(QScrollArea)
            if scroll:
                scroll.ensureWidgetVisible(targets[idx], 0, 20)
        for i, btn in enumerate(self.nav_buttons):
            btn.setProperty("active", "true" if i == idx else "false")
            btn.style().unpolish(btn)
            btn.style().polish(btn)

    def _add_panel(self, groupbox):
        self.panels_layout.addWidget(groupbox)

    # -----------------------------------------------------------------
    # Panel 1: Targets & Output
    # -----------------------------------------------------------------

    def _create_panel_targets_output(self):
        self.grp_targets = QGroupBox("1. Targets && Output")
        lay = QVBoxLayout(self.grp_targets)
        lay.setSpacing(8)

        self.pick_targets_dir = FilePickerRow(
            "Targets Directory:", placeholder="Folder containing the FITS datacubes",
            is_directory=True
        )
        lay.addWidget(self.pick_targets_dir)

        self.pick_targets_csv = FilePickerRow(
            "Targets CSV File:", placeholder="Path where the table below is written before each run",
            file_filter="CSV Files (*.csv);;All Files (*)"
        )
        lay.addWidget(self.pick_targets_csv)

        self.targets_table = TargetsTableWidget()
        self.targets_table.csv_path_changed.connect(self.pick_targets_csv.set_text)
        lay.addWidget(self.targets_table)

        self.pick_output_dir = FilePickerRow(
            "Output Directory:", placeholder="./urutau_output/", is_directory=True
        )
        self.pick_output_dir.set_text("./urutau_output/")
        lay.addWidget(self.pick_output_dir)

        row = QWidget()
        row_lay = QHBoxLayout(row)
        row_lay.setContentsMargins(0, 0, 0, 0)
        row_lay.setSpacing(16)
        self.chk_save_config = QCheckBox("  Save config HDU with each result")
        self.chk_save_config.setChecked(True)
        self.chk_debug = QCheckBox("  Debug output (print module parameters)")
        self.chk_debug.setChecked(True)
        row_lay.addWidget(self.chk_save_config)
        row_lay.addWidget(self.chk_debug)
        row_lay.addStretch()
        lay.addWidget(row)

        self.chk_overwrite = QCheckBox("  Overwrite existing output files")
        self.chk_overwrite.setChecked(True)
        lay.addWidget(self.chk_overwrite)

        info_overwrite = QLabel(
            "Unchecked: a target whose megacube already exists in the Output Directory is "
            "skipped entirely (none of its modules run, including STARLIGHT) instead of "
            "being recomputed and overwritten."
        )
        info_overwrite.setProperty("muted", "true")
        info_overwrite.setWordWrap(True)
        lay.addWidget(info_overwrite)

        threads_row, self.spin_num_threads = labeled_int(
            "Urutau Threads:", min_val=1, max_val=256, default_val=1,
            suffix="parallel targets"
        )
        lay.addWidget(threads_row)

        self._add_panel(self.grp_targets)

    # -----------------------------------------------------------------
    # Panel 2: Input HDUs
    # -----------------------------------------------------------------

    def _create_panel_input_hdus(self):
        self.grp_input = QGroupBox("2. Input HDUs")
        lay = QVBoxLayout(self.grp_input)
        lay.setSpacing(8)

        row, self.edit_data_hdu = labeled_text("Data (flux) HDU:", default_val="DATA")
        lay.addWidget(row)
        row, self.edit_stat_hdu = labeled_text("Stat/Error HDU:", default_val="STAT")
        lay.addWidget(row)
        row, self.combo_stat_type = labeled_combo(
            "Stat HDU Type:", list(STAT_TYPE_MAP.keys()), default_index=0
        )
        lay.addWidget(row)

        info = QLabel(
            "These are the starting extensions in the input FITS files. Every "
            "enabled stage below is automatically chained onto its predecessor's "
            "output HDU (matching Urutau's *_RSP / *_DEGR / *_BIN naming)."
        )
        info.setProperty("muted", "true")
        info.setWordWrap(True)
        lay.addWidget(info)

        self._add_panel(self.grp_input)

    # -----------------------------------------------------------------
    # Panel 3: Spatial Resampling
    # -----------------------------------------------------------------

    def _create_panel_spatial(self):
        self.grp_spatial = QGroupBox("3. Spatial Resampling  (optional)")
        self.grp_spatial.setCheckable(True)
        self.grp_spatial.setChecked(False)
        outer = QVBoxLayout(self.grp_spatial)
        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        row, self.spin_spatial_size = labeled_int(
            "Bin Size:", min_val=1, max_val=1000, default_val=4, suffix="pixels/side"
        )
        lay.addWidget(row)

        info = QLabel("Bins the datacube spatially into NxN spaxel blocks before spectral processing.")
        info.setProperty("muted", "true")
        info.setWordWrap(True)
        lay.addWidget(info)

        outer.addWidget(content)
        CheckableGroupBox.wire(self.grp_spatial, content)
        self._add_panel(self.grp_spatial)

    # -----------------------------------------------------------------
    # Panel 4: Resolution degradation
    # -----------------------------------------------------------------

    def _create_panel_resolution(self):
        self.grp_resolution = QGroupBox("4. Spectral Resolution Degradation  (optional)")
        self.grp_resolution.setCheckable(True)
        self.grp_resolution.setChecked(True)
        outer = QVBoxLayout(self.grp_resolution)
        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        units = ["R", "FWHM", "sigma"]

        row, self.combo_res_input_type = labeled_combo("Input Unit:", units, default_index=0)
        lay.addWidget(row)
        row, self.spin_res_input_value = labeled_double(
            "Input Value:", min_val=0., max_val=1e7, decimals=3, default_val=3027.0
        )
        lay.addWidget(row)
        self.opt_res_input_ref = OptionalDoubleRow(
            "Input Ref. λ (optional)", suffix="Å", max_val=1e6, decimals=2, default_val=5500.0
        )
        lay.addWidget(self.opt_res_input_ref)

        row, self.combo_res_output_type = labeled_combo("Output Unit:", units, default_index=1)
        lay.addWidget(row)
        row, self.spin_res_output_value = labeled_double(
            "Output Value:", min_val=0., max_val=1e7, decimals=4, default_val=2.51
        )
        lay.addWidget(row)
        self.opt_res_output_ref = OptionalDoubleRow(
            "Output Ref. λ (optional)", suffix="Å", max_val=1e6, decimals=2, default_val=5500.0
        )
        lay.addWidget(self.opt_res_output_ref)

        info = QLabel(
            "Applied to both the data and stat HDUs (DegradeDataFlex). Convolution can only "
            "broaden resolution: the output must be coarser than the input."
        )
        info.setProperty("muted", "true")
        info.setWordWrap(True)
        lay.addWidget(info)

        outer.addWidget(content)
        CheckableGroupBox.wire(self.grp_resolution, content)
        self._add_panel(self.grp_resolution)

    # -----------------------------------------------------------------
    # Panel 5: Spectral binning
    # -----------------------------------------------------------------

    def _create_panel_spectral_bin(self):
        self.grp_spectral_bin = QGroupBox("5. Spectral Resampling / Binning  (optional)")
        self.grp_spectral_bin.setCheckable(True)
        self.grp_spectral_bin.setChecked(True)
        outer = QVBoxLayout(self.grp_spectral_bin)
        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        row, self.spin_spectral_bin_size = labeled_double(
            "Sample Size:", min_val=0.0001, max_val=1000.0, decimals=4, default_val=1.0,
            suffix="Å/pixel"
        )
        lay.addWidget(row)

        outer.addWidget(content)
        CheckableGroupBox.wire(self.grp_spectral_bin, content)
        self._add_panel(self.grp_spectral_bin)

    # -----------------------------------------------------------------
    # Panel 6: Dereddening
    # -----------------------------------------------------------------

    def _create_panel_dereddening(self):
        self.grp_dereddening = QGroupBox("6. Galactic Dereddening  (optional)")
        self.grp_dereddening.setCheckable(True)
        self.grp_dereddening.setChecked(True)
        outer = QVBoxLayout(self.grp_dereddening)
        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        row, self.combo_dered_law = labeled_combo(
            "Law:", list(DEREDDENING_LAWS.keys()), default_index=0
        )
        lay.addWidget(row)
        row, self.spin_dered_ebv = labeled_double(
            "E(B-V):", min_val=0.0, max_val=10.0, decimals=4, default_val=0.0
        )
        lay.addWidget(row)
        row, self.spin_dered_rv = labeled_double(
            "Rv:", min_val=0.0, max_val=10.0, decimals=3, default_val=3.1
        )
        lay.addWidget(row)

        info = QLabel("Output flux HDU is always named FLUX_DRD (fixed by Urutau).")
        info.setProperty("muted", "true")
        lay.addWidget(info)

        outer.addWidget(content)
        CheckableGroupBox.wire(self.grp_dereddening, content)
        self._add_panel(self.grp_dereddening)

    # -----------------------------------------------------------------
    # Panel 7: S/N Mask
    # -----------------------------------------------------------------

    def _create_panel_sn_mask(self):
        self.grp_sn_mask = QGroupBox("7. Signal-to-Noise Mask  (optional, feeds Starlight's flag HDU)")
        self.grp_sn_mask.setCheckable(True)
        self.grp_sn_mask.setChecked(True)
        outer = QVBoxLayout(self.grp_sn_mask)
        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        row, self.combo_sn_method = labeled_combo(
            "Method:", list(SN_MASK_METHODS.keys()), default_index=0
        )
        lay.addWidget(row)

        window_row = QWidget()
        window_lay = QHBoxLayout(window_row)
        window_lay.setContentsMargins(0, 0, 0, 0)
        window_lay.setSpacing(8)
        window_label = QLabel("S/N Window (Å):")
        window_label.setFixedWidth(150)
        window_label.setStyleSheet("font-weight: 500;")
        window_lay.addWidget(window_label)
        self.spin_sn_window_min = self._plain_double(0., 1e6, 2, 5650.0)
        self.spin_sn_window_max = self._plain_double(0., 1e6, 2, 5750.0)
        window_lay.addWidget(self.spin_sn_window_min, 1)
        window_lay.addWidget(QLabel("to"))
        window_lay.addWidget(self.spin_sn_window_max, 1)
        lay.addWidget(window_row)

        row, self.edit_sn_thresholds = labeled_text(
            "Thresholds:", placeholder="Comma-separated, e.g. 10,20", default_val="10,20"
        )
        lay.addWidget(row)

        row, self.spin_sn_redshift = labeled_double(
            "Redshift (window):", min_val=0.0, max_val=10.0, decimals=6, default_val=0.0
        )
        lay.addWidget(row)

        outer.addWidget(content)
        CheckableGroupBox.wire(self.grp_sn_mask, content)
        self._add_panel(self.grp_sn_mask)

    def _plain_double(self, min_val, max_val, decimals, default_val):
        from PyQt5.QtWidgets import QDoubleSpinBox
        spin = QDoubleSpinBox()
        spin.setRange(min_val, max_val)
        spin.setDecimals(decimals)
        spin.setValue(default_val)
        return spin

    # -----------------------------------------------------------------
    # Panel 8: Starlight
    # -----------------------------------------------------------------

    def _create_panel_starlight(self):
        self.grp_starlight = QGroupBox("8. Starlight — Stellar Population Fitting  (optional)")
        self.grp_starlight.setCheckable(True)
        self.grp_starlight.setChecked(True)
        outer = QVBoxLayout(self.grp_starlight)
        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)

        self.pick_starlight_exe = FilePickerRow(
            "Executable Path:", placeholder="./starlight/StarlightChains_v04...static.exe"
        )
        lay.addWidget(self.pick_starlight_exe)

        self.pick_starlight_grid = FilePickerRow(
            "Default Grid File:", placeholder="./starlight/reference_grid.in"
        )
        lay.addWidget(self.pick_starlight_grid)

        row, self.edit_mask_file = labeled_text(
            "Mask File Override:", placeholder="Optional — leave empty to use the grid file's own mask"
        )
        lay.addWidget(row)

        info_mask = QLabel(
            "Without an override, the mask filename comes from the Default Grid File above. "
            "The field here replaces it for every target. To give each target its own "
            "emission-line mask instead (or on top of the field above), add a 'mask file' "
            "column to the Targets table (Section 1) with each mask's filename (must live "
            "in the grid's mask directory) — a per-target CSV value always wins over this field."
        )
        info_mask.setProperty("muted", "true")
        info_mask.setWordWrap(True)
        lay.addWidget(info_mask)

        row, self.spin_starlight_threads = labeled_int(
            "STARLIGHT Threads:", min_val=1, max_val=512, default_val=1,
            suffix="parallel executables"
        )
        lay.addWidget(row)

        row, self.spin_flag_threshold = labeled_int(
            "S/N Flag Threshold:", min_val=0, max_val=100000, default_val=10,
            suffix="(must match an S/N Mask threshold)"
        )
        lay.addWidget(row)

        row, self.spin_gal_distance = labeled_double(
            "Galaxy Distance:", min_val=0.0, max_val=1e6, decimals=3, default_val=0.0,
            suffix="Mpc"
        )
        lay.addWidget(row)
        row, self.spin_starlight_redshift = labeled_double(
            "Redshift:", min_val=0.0, max_val=10.0, decimals=6, default_val=0.0
        )
        lay.addWidget(row)
        row, self.spin_norm_factor = labeled_double(
            "Normalization Factor:", min_val=1e-30, max_val=1e30, decimals=10, default_val=1.0
        )
        lay.addWidget(row)
        row, self.edit_flux_unit = labeled_text("Flux Unit:", placeholder="e.g. erg/s/cm2/A")
        lay.addWidget(row)

        self.chk_keep_tmp = QCheckBox("  Keep STARLIGHT temporary files")
        lay.addWidget(self.chk_keep_tmp)

        lay.addWidget(self._section_label("Per-Spaxel Process Timeout  (optional)"))

        row, self.combo_timeout_mode = labeled_combo(
            "Timeout Mode:",
            ["None (never kill)", "Fixed", "Adaptive (rolling average)"],
            default_index=0,
        )
        self.combo_timeout_mode.currentIndexChanged.connect(self._on_timeout_mode_changed)
        lay.addWidget(row)

        row, self.spin_timeout_minutes = labeled_double(
            "Timeout Minutes:", min_val=0.0, max_val=10000.0, decimals=2, default_val=15.0,
            suffix="min — 0 = no limit"
        )
        lay.addWidget(row)

        row, self.spin_timeout_window = labeled_int(
            "Rolling Window:", min_val=2, max_val=10000, default_val=15,
            suffix="recent spaxels averaged"
        )
        lay.addWidget(row)

        row, self.spin_timeout_multiplier = labeled_double(
            "Timeout Multiplier:", min_val=1.0, max_val=100.0, decimals=2, default_val=2.0,
            suffix="× the rolling average"
        )
        lay.addWidget(row)

        info_timeout = QLabel(
            "Fixed: kills any spaxel's STARLIGHT process past 'Timeout Minutes'. Adaptive: "
            "kills a process running longer than 'Timeout Multiplier' × the average duration "
            "of the last 'Rolling Window' spaxels that finished normally in this same target "
            "— 'Timeout Minutes' is used as a fallback limit until enough of them have run. A "
            "killed spaxel is simply recorded as a failed spaxel, like any other STARLIGHT "
            "failure."
        )
        info_timeout.setProperty("muted", "true")
        info_timeout.setWordWrap(True)
        lay.addWidget(info_timeout)

        self._on_timeout_mode_changed()

        pop_header = QWidget()
        pop_header_lay = QHBoxLayout(pop_header)
        pop_header_lay.setContentsMargins(0, 0, 0, 0)
        pop_header_lay.addWidget(self._section_label(
            "Population Bins  (required — name: [min, max] age in years, "
            "min < age ≤ max)"
        ))
        pop_header_lay.addStretch()
        btn_check_agn = QPushButton("🔍 Check vs Base (AGN conflicts)")
        btn_check_agn.setObjectName("minorBtn")
        btn_check_agn.setCursor(Qt.PointingHandCursor)
        btn_check_agn.clicked.connect(self._on_check_base_agn)
        pop_header_lay.addWidget(btn_check_agn)
        lay.addWidget(pop_header)

        info_agn = QLabel(
            "Cross-checks these bins against the base file referenced by the default grid "
            "file (resolved next to the Starlight executable). Any base component whose name "
            "hints at AGN must fall outside every bin below — it belongs in the "
            "Featureless-Continuum or Black-Body tables instead."
        )
        info_agn.setProperty("muted", "true")
        info_agn.setWordWrap(True)
        lay.addWidget(info_agn)

        self.table_population_ages = KeyRangeTable(min_header="Min (excl.)", max_header="Max (incl.)")
        lay.addWidget(self.table_population_ages)

        lay.addWidget(self._section_label(
            "SFR Ages  (optional — name: [min, max] age in years, min < age ≤ max)"
        ))
        self.table_sfr_ages = KeyRangeTable(min_header="Min (excl.)", max_header="Max (incl.)")
        lay.addWidget(self.table_sfr_ages)

        lay.addWidget(self._section_label(
            "Returned-Mass Ages  (optional — name: [min, max] age in years, min < age ≤ max)"
        ))
        self.table_ret_mass_ages = KeyRangeTable(min_header="Min (excl.)", max_header="Max (incl.)")
        lay.addWidget(self.table_ret_mass_ages)

        lay.addWidget(self._section_label(
            "Featureless-Continuum (AGN) Exponents  (optional — name: [min, max] exponent, "
            "min < exp ≤ max)"
        ))
        self.table_fc_exps = KeyRangeTable(min_header="Min (excl.)", max_header="Max (incl.)")
        lay.addWidget(self.table_fc_exps)

        lay.addWidget(self._section_label(
            "Black-Body Temperatures  (optional — name: [min, max] temperature in K, "
            "min < temp ≤ max)"
        ))
        self.table_bb_temps = KeyRangeTable(min_header="Min (excl.)", max_header="Max (incl.)")
        lay.addWidget(self.table_bb_temps)

        outer.addWidget(content)
        CheckableGroupBox.wire(self.grp_starlight, content)
        self._add_panel(self.grp_starlight)

    def _section_label(self, text):
        lbl = QLabel(text)
        lbl.setStyleSheet(f"color: {ACCENT}; font-weight: 600; font-size: 12px; margin-top: 6px;")
        return lbl

    def _on_timeout_mode_changed(self):
        mode = self._timeout_mode_key()
        self.spin_timeout_minutes.setEnabled(mode in ("fixed", "adaptive"))
        self.spin_timeout_window.setEnabled(mode == "adaptive")
        self.spin_timeout_multiplier.setEnabled(mode == "adaptive")

    def _timeout_mode_key(self) -> str:
        return {
            "None (never kill)": "none",
            "Fixed": "fixed",
            "Adaptive (rolling average)": "adaptive",
        }[self.combo_timeout_mode.currentText()]

    def _on_check_base_agn(self):
        grid_file = self.pick_starlight_grid.text()
        starlight_path = self.pick_starlight_exe.text()

        try:
            components = read_base_components(grid_file, starlight_path)
        except BaseGridError as e:
            QMessageBox.warning(self, "Can't Read Base File", str(e))
            return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Unexpected error reading base file: {e}")
            return

        try:
            population_bins = self.table_population_ages.get_dict()
        except ValueError as e:
            QMessageBox.critical(self, "Invalid Population Bins", str(e))
            return

        agn_components = [c for c in components if c["is_agn"]]
        conflicts = find_agn_bin_conflicts(components, population_bins)

        lines = [f"Base has {len(components)} components ({len(agn_components)} AGN-hinted)."]
        if agn_components:
            lines.append("")
            lines.append("AGN-hinted components:")
            for c in agn_components:
                lines.append(f"  • {c['name']}  (age = {c['age']:.4g} yr)")

        if conflicts:
            lines.append("")
            lines.append("⚠ Remove these from the Population Bins table — they belong in the "
                          "Featureless-Continuum / Black-Body tables instead:")
            for conf in conflicts:
                lines.append(
                    f"  • '{conf['component']}' (age = {conf['age']:.4g} yr) falls inside "
                    f"bin '{conf['bin']}'"
                )
            QMessageBox.warning(self, "AGN / Population Bin Conflicts", "\n".join(lines))
        else:
            lines.append("")
            lines.append("No AGN components fall inside the current Population Bins.")
            QMessageBox.information(self, "Base Cross-Check", "\n".join(lines))

    # -----------------------------------------------------------------
    # Bottom Panel: Progress + Log
    # -----------------------------------------------------------------

    def _create_bottom_panel(self):
        bottom = QWidget()
        lay = QVBoxLayout(bottom)
        lay.setContentsMargins(20, 8, 20, 12)
        lay.setSpacing(8)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        lay.addWidget(self.progress)

        lbl_log = QLabel("Run Log")
        lbl_log.setStyleSheet(f"color: {ACCENT}; font-size: 12px; font-weight: bold;")
        lay.addWidget(lbl_log)

        self.log_console = LogConsole()
        lay.addWidget(self.log_console, 1)

        return bottom

    # -----------------------------------------------------------------
    # Config gather / validate / apply / run
    # -----------------------------------------------------------------

    def _gather_config(self) -> dict:
        cfg = {
            "num_threads": self.spin_num_threads.value(),
            "input": {
                "data_hdu": self.edit_data_hdu.text().strip() or "DATA",
                "stat_hdu": self.edit_stat_hdu.text().strip() or "STAT",
                "stat_type": self.combo_stat_type.currentText(),
            },
            "spatial": {
                "enabled": self.grp_spatial.isChecked(),
                "size": self.spin_spatial_size.value(),
            },
            "resolution": {
                "enabled": self.grp_resolution.isChecked(),
                "input_type": self.combo_res_input_type.currentText(),
                "input_value": self.spin_res_input_value.value(),
                "input_ref_wave": self.opt_res_input_ref.value(),
                "output_type": self.combo_res_output_type.currentText(),
                "output_value": self.spin_res_output_value.value(),
                "output_ref_wave": self.opt_res_output_ref.value(),
            },
            "spectral_bin": {
                "enabled": self.grp_spectral_bin.isChecked(),
                "size": self.spin_spectral_bin_size.value(),
            },
            "dereddening": {
                "enabled": self.grp_dereddening.isChecked(),
                "law": self.combo_dered_law.currentText(),
                "ebv": self.spin_dered_ebv.value(),
                "rv": self.spin_dered_rv.value(),
            },
            "sn_mask": {
                "enabled": self.grp_sn_mask.isChecked(),
                "method": self.combo_sn_method.currentText(),
                "window": (self.spin_sn_window_min.value(), self.spin_sn_window_max.value()),
                "thresholds": self._parse_thresholds(self.edit_sn_thresholds.text()),
                "redshift": self.spin_sn_redshift.value(),
            },
            "starlight": {
                "enabled": self.grp_starlight.isChecked(),
                "path": self.pick_starlight_exe.text(),
                "grid_file": self.pick_starlight_grid.text(),
                "num_threads": self.spin_starlight_threads.value(),
                "flag_threshold": self.spin_flag_threshold.value(),
                "galaxy_distance": self.spin_gal_distance.value(),
                "redshift": self.spin_starlight_redshift.value(),
                "normalization_factor": self.spin_norm_factor.value(),
                "flux_unit": self.edit_flux_unit.text().strip(),
                "keep_tmp": self.chk_keep_tmp.isChecked(),
                "mask_file": self.edit_mask_file.text().strip(),
                "timeout_mode": self._timeout_mode_key(),
                "timeout_minutes": (
                    self.spin_timeout_minutes.value()
                    if self._timeout_mode_key() != "none" and self.spin_timeout_minutes.value() > 0
                    else None
                ),
                "timeout_window": self.spin_timeout_window.value(),
                "timeout_multiplier": self.spin_timeout_multiplier.value(),
                "population_ages": self.table_population_ages.get_dict(),
                "sfr_ages": self.table_sfr_ages.get_dict(),
                "ret_mass_ages": self.table_ret_mass_ages.get_dict(),
                "fc_exps": self.table_fc_exps.get_dict(),
                "bb_temps": self.table_bb_temps.get_dict(),
            },
            "targets": self._gather_targets_config(),
            "output": {
                "save_path_root": self.pick_output_dir.text() or "./urutau_output/",
                "save_config": self.chk_save_config.isChecked(),
                "debug": self.chk_debug.isChecked(),
                "overwrite": self.chk_overwrite.isChecked(),
            },
        }
        return cfg

    def _gather_targets_config(self) -> dict:
        columns, rows = self.targets_table.get_columns_rows()
        return {
            "dir": self.pick_targets_dir.text(),
            "csv": self.pick_targets_csv.text(),
            "columns": columns,
            "rows": rows,
        }

    @staticmethod
    def _parse_thresholds(text):
        values = []
        for tok in text.split(","):
            tok = tok.strip()
            if not tok:
                continue
            values.append(float(tok) if "." in tok else int(tok))
        return values

    def _validate(self, cfg, check_paths=True) -> list:
        """
        check_paths=False skips local-filesystem existence checks (Targets
        Directory / CSV) — used when exporting a standalone script that may
        be intended to run later, or on a different machine/mount.
        """
        errors = []

        targets = cfg["targets"]
        if not targets["dir"]:
            errors.append("Targets Directory is required.")
        elif check_paths and not os.path.isdir(targets["dir"]):
            errors.append("Targets Directory must point to an existing folder.")
        if not targets["csv"]:
            errors.append("Targets CSV File needs a path (the table is written there before running).")
        if not any(row and str(row[0]).strip() for row in targets["rows"]):
            errors.append("Targets table needs at least one row with a target filename in the first column.")

        if not any([
            cfg["spatial"]["enabled"], cfg["resolution"]["enabled"],
            cfg["spectral_bin"]["enabled"], cfg["dereddening"]["enabled"],
            cfg["sn_mask"]["enabled"], cfg["starlight"]["enabled"],
        ]):
            errors.append("Enable at least one processing stage.")

        if cfg["sn_mask"]["enabled"] and not cfg["sn_mask"]["thresholds"]:
            errors.append("S/N Mask: provide at least one threshold.")

        if cfg["starlight"]["enabled"]:
            sl = cfg["starlight"]
            if not sl["path"]:
                errors.append("Starlight: executable path is required.")
            if not sl["grid_file"]:
                errors.append("Starlight: default grid file is required.")
            if not sl["population_ages"]:
                errors.append("Starlight: at least one Population Bin row is required.")
            if cfg["sn_mask"]["enabled"] and sl["flag_threshold"] not in cfg["sn_mask"]["thresholds"]:
                errors.append(
                    "Starlight: S/N Flag Threshold must match one of the S/N Mask thresholds."
                )
            if sl["timeout_mode"] == "fixed" and not sl["timeout_minutes"]:
                errors.append(
                    "Starlight: Timeout Mode is 'Fixed' but Timeout Minutes is 0 — "
                    "set a value greater than 0, or switch the mode to 'None'."
                )

            # Best-effort cross-check against the base file: only enforced when the
            # grid/base can actually be resolved and read (e.g. not on a remote mount).
            if sl["path"] and sl["grid_file"] and sl["population_ages"]:
                try:
                    components = read_base_components(sl["grid_file"], sl["path"])
                    conflicts = find_agn_bin_conflicts(components, sl["population_ages"])
                except BaseGridError:
                    conflicts = None
                except Exception:
                    conflicts = None
                if conflicts:
                    for conf in conflicts:
                        errors.append(
                            f"Starlight: AGN component '{conf['component']}' (age = "
                            f"{conf['age']:.4g} yr) falls inside Population Bin '{conf['bin']}' — "
                            "move it to the Featureless-Continuum or Black-Body table instead."
                        )

        return errors

    def _on_run(self):
        if self.worker and self.worker.isRunning():
            QMessageBox.warning(self, "Running", "Urutau is already running.")
            return

        try:
            cfg = self._gather_config()
        except ValueError as e:
            QMessageBox.critical(self, "Invalid Parameter", str(e))
            return

        errors = self._validate(cfg)
        if errors:
            QMessageBox.warning(self, "Missing / Invalid Input", "\n".join(f"• {e}" for e in errors))
            return

        try:
            self.targets_table.save_csv(cfg["targets"]["csv"])
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not write targets CSV: {e}")
            return

        self.log_console.clear_log()
        self.progress.setVisible(True)
        self.btn_run.setEnabled(False)
        self.btn_run.setText("⏳  Running…")
        self.statusBar().showMessage("Running Urutau…")

        self.worker = UrutauWorker(cfg)
        self.worker.log_signal.connect(self.log_console.append_log)
        self.worker.finished_signal.connect(self._on_worker_finished)
        self.worker.start()

    def _on_worker_finished(self, success):
        self.progress.setVisible(False)
        self.btn_run.setEnabled(True)
        self.btn_run.setText("▶  Run Urutau")
        if success:
            self.log_console.append_log("\n✅ Urutau finished successfully.")
            self.statusBar().showMessage("Done.")
        else:
            self.statusBar().showMessage("Run finished with errors. Check the log.")

    # -----------------------------------------------------------------
    # Save / Load Config
    # -----------------------------------------------------------------

    def _apply_config(self, cfg: dict):
        self.spin_num_threads.setValue(cfg.get("num_threads", 1))

        inp = cfg.get("input", {})
        self.edit_data_hdu.setText(inp.get("data_hdu", "DATA"))
        self.edit_stat_hdu.setText(inp.get("stat_hdu", "STAT"))
        idx = self.combo_stat_type.findText(inp.get("stat_type", "variance"))
        self.combo_stat_type.setCurrentIndex(max(idx, 0))

        sp = cfg.get("spatial", {})
        self.grp_spatial.setChecked(sp.get("enabled", False))
        self.spin_spatial_size.setValue(sp.get("size", 4))

        res = cfg.get("resolution", {})
        self.grp_resolution.setChecked(res.get("enabled", True))
        self.combo_res_input_type.setCurrentText(res.get("input_type", "R"))
        self.spin_res_input_value.setValue(res.get("input_value", 3027.0))
        self.opt_res_input_ref.set_state(
            res.get("input_ref_wave") is not None, res.get("input_ref_wave")
        )
        self.combo_res_output_type.setCurrentText(res.get("output_type", "FWHM"))
        self.spin_res_output_value.setValue(res.get("output_value", 2.51))
        self.opt_res_output_ref.set_state(
            res.get("output_ref_wave") is not None, res.get("output_ref_wave")
        )

        sb = cfg.get("spectral_bin", {})
        self.grp_spectral_bin.setChecked(sb.get("enabled", True))
        self.spin_spectral_bin_size.setValue(sb.get("size", 1.0))

        dr = cfg.get("dereddening", {})
        self.grp_dereddening.setChecked(dr.get("enabled", True))
        self.combo_dered_law.setCurrentText(dr.get("law", "CCM"))
        self.spin_dered_ebv.setValue(dr.get("ebv", 0.0))
        self.spin_dered_rv.setValue(dr.get("rv", 3.1))

        sn = cfg.get("sn_mask", {})
        self.grp_sn_mask.setChecked(sn.get("enabled", True))
        self.combo_sn_method.setCurrentText(sn.get("method", "Variance"))
        window = sn.get("window", (5650.0, 5750.0))
        self.spin_sn_window_min.setValue(window[0])
        self.spin_sn_window_max.setValue(window[1])
        self.edit_sn_thresholds.setText(
            ",".join(str(t) for t in sn.get("thresholds", [10, 20]))
        )
        self.spin_sn_redshift.setValue(sn.get("redshift", 0.0))

        sl = cfg.get("starlight", {})
        self.grp_starlight.setChecked(sl.get("enabled", True))
        self.pick_starlight_exe.set_text(sl.get("path", ""))
        self.pick_starlight_grid.set_text(sl.get("grid_file", ""))
        self.spin_starlight_threads.setValue(sl.get("num_threads", 1))
        self.spin_flag_threshold.setValue(sl.get("flag_threshold") or 0)
        self.spin_gal_distance.setValue(sl.get("galaxy_distance", 0.0))
        self.spin_starlight_redshift.setValue(sl.get("redshift", 0.0))
        self.spin_norm_factor.setValue(sl.get("normalization_factor", 1.0))
        self.edit_flux_unit.setText(sl.get("flux_unit", ""))
        self.chk_keep_tmp.setChecked(sl.get("keep_tmp", False))
        self.edit_mask_file.setText(sl.get("mask_file", ""))
        mode_labels = {
            "none": "None (never kill)",
            "fixed": "Fixed",
            "adaptive": "Adaptive (rolling average)",
        }
        self.combo_timeout_mode.setCurrentText(mode_labels.get(sl.get("timeout_mode", "none"), "None (never kill)"))
        self.spin_timeout_minutes.setValue(sl.get("timeout_minutes") or 15.0)
        self.spin_timeout_window.setValue(sl.get("timeout_window", 15))
        self.spin_timeout_multiplier.setValue(sl.get("timeout_multiplier", 2.0))
        self._on_timeout_mode_changed()
        self.table_population_ages.set_dict(sl.get("population_ages", {}))
        self.table_sfr_ages.set_dict(sl.get("sfr_ages", {}))
        self.table_ret_mass_ages.set_dict(sl.get("ret_mass_ages", {}))
        self.table_fc_exps.set_dict(sl.get("fc_exps", {}))
        self.table_bb_temps.set_dict(sl.get("bb_temps", {}))

        tg = cfg.get("targets", {})
        self.pick_targets_dir.set_text(tg.get("dir", ""))
        self.pick_targets_csv.set_text(tg.get("csv", ""))
        self.targets_table.set_columns_rows(tg.get("columns", ["target"]), tg.get("rows", []))

        out = cfg.get("output", {})
        self.pick_output_dir.set_text(out.get("save_path_root", "./urutau_output/"))
        self.chk_save_config.setChecked(out.get("save_config", True))
        self.chk_debug.setChecked(out.get("debug", True))
        self.chk_overwrite.setChecked(out.get("overwrite", True))

    def _on_save_config(self):
        try:
            cfg = self._gather_config()
        except ValueError as e:
            QMessageBox.critical(self, "Invalid Parameter", str(e))
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Urutau Configuration", "urutau_config.json",
            "JSON Files (*.json);;All Files (*)"
        )
        if path:
            with open(path, "w") as f:
                json.dump(cfg, f, indent=2, default=str)
            self.statusBar().showMessage(f"Configuration saved to {path}")

    def _on_load_config(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Urutau Configuration", "", "JSON Files (*.json);;All Files (*)"
        )
        if path:
            with open(path) as f:
                cfg = json.load(f)
            self._apply_config(cfg)
            self.statusBar().showMessage(f"Configuration loaded from {path}")

    def _on_export_script(self):
        """
        Writes the current configuration out as a standalone .py script that
        reproduces the same Urutau pipeline without needing the GUI — meant
        to be run by hand (e.g. on a cluster with no display).
        """
        try:
            cfg = self._gather_config()
        except ValueError as e:
            QMessageBox.critical(self, "Invalid Parameter", str(e))
            return

        # Local-path existence isn't required: the script may target a
        # different machine/mount than the one the GUI is running on.
        errors = self._validate(cfg, check_paths=False)
        if errors:
            QMessageBox.warning(self, "Missing / Invalid Input", "\n".join(f"• {e}" for e in errors))
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Urutau Script", "run_urutau.py", "Python Files (*.py);;All Files (*)"
        )
        if not path:
            return

        if cfg["targets"]["csv"]:
            try:
                self.targets_table.save_csv(cfg["targets"]["csv"])
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not write targets CSV: {e}")
                return

        try:
            script = generate_script(cfg)
            with open(path, "w") as f:
                f.write(script)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not write script: {e}")
            return

        self.statusBar().showMessage(f"Script exported to {path}")
        QMessageBox.information(
            self, "Script Exported",
            f"Saved to:\n{path}\n\nThe targets CSV was (re)written to:\n{cfg['targets']['csv']}\n\n"
            f"Run it with:\n    python {os.path.basename(path)}"
        )


# ---------------------------------------------------------------------------
# Entry Point
# ---------------------------------------------------------------------------

def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setApplicationName("Urutau")

    font = QFont("Inter", 10)
    app.setFont(font)

    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
