"""
custom_widgets.py — Reusable GUI widgets for the Urutau interface.
"""

import os

import pandas as pd

from PyQt5.QtWidgets import (
    QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QLabel, QLineEdit,
    QPushButton, QFileDialog, QCheckBox, QDoubleSpinBox, QComboBox,
    QSpinBox, QTextEdit, QSizePolicy, QTableWidget, QTableWidgetItem,
    QHeaderView, QInputDialog, QMessageBox
)
from PyQt5.QtCore import Qt, pyqtSignal


LABEL_WIDTH = 150


class FilePickerRow(QWidget):
    """Label + QLineEdit + Browse button, for files or directories."""

    def __init__(self, label_text, placeholder="", is_directory=False,
                 file_filter="All Files (*)", parent=None):
        super().__init__(parent)
        self.is_directory = is_directory
        self.file_filter = file_filter

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self.label = QLabel(label_text)
        self.label.setFixedWidth(LABEL_WIDTH)
        self.label.setStyleSheet("font-weight: 500;")
        layout.addWidget(self.label)

        self.line_edit = QLineEdit()
        self.line_edit.setPlaceholderText(placeholder)
        self.line_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        layout.addWidget(self.line_edit, 1)

        self.btn_browse = QPushButton("Browse…")
        self.btn_browse.setObjectName("minorBtn")
        self.btn_browse.setFixedWidth(90)
        self.btn_browse.setCursor(Qt.PointingHandCursor)
        self.btn_browse.clicked.connect(self._on_browse)
        layout.addWidget(self.btn_browse)

    def _on_browse(self):
        if self.is_directory:
            path = QFileDialog.getExistingDirectory(
                self, f"Select {self.label.text().strip(':')}", self.line_edit.text()
            )
        else:
            path, _ = QFileDialog.getOpenFileName(
                self, f"Select {self.label.text().strip(':')}", self.line_edit.text(),
                self.file_filter
            )
        if path:
            self.line_edit.setText(path)

    def text(self):
        return self.line_edit.text().strip()

    def set_text(self, value):
        self.line_edit.setText(value or "")


class LabeledRow(QWidget):
    """Generic label + single child widget row."""

    def __init__(self, label_text, child, label_width=LABEL_WIDTH, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        label = QLabel(label_text)
        label.setFixedWidth(label_width)
        label.setStyleSheet("font-weight: 500;")
        layout.addWidget(label)

        child.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        layout.addWidget(child, 1)
        self.child = child


def labeled_double(label_text, min_val=-1e9, max_val=1e9, decimals=4,
                   default_val=0.0, suffix=""):
    spin = QDoubleSpinBox()
    spin.setRange(min_val, max_val)
    spin.setDecimals(decimals)
    spin.setValue(default_val)
    if suffix:
        spin.setSuffix(f"  {suffix}")
    row = LabeledRow(label_text, spin)
    return row, spin


def labeled_int(label_text, min_val=0, max_val=10_000_000, default_val=1, suffix=""):
    spin = QSpinBox()
    spin.setRange(min_val, max_val)
    spin.setValue(default_val)
    if suffix:
        spin.setSuffix(f"  {suffix}")
    row = LabeledRow(label_text, spin)
    return row, spin


def labeled_combo(label_text, options, default_index=0):
    combo = QComboBox()
    combo.addItems(options)
    combo.setCurrentIndex(default_index)
    row = LabeledRow(label_text, combo)
    return row, combo


def labeled_text(label_text, placeholder="", default_val=""):
    edit = QLineEdit()
    edit.setPlaceholderText(placeholder)
    edit.setText(default_val)
    row = LabeledRow(label_text, edit)
    return row, edit


class OptionalDoubleRow(QWidget):
    """A checkbox-controlled optional double value (e.g. an optional reference wavelength)."""

    def __init__(self, label_text, suffix="", min_val=-1e9, max_val=1e9,
                 decimals=2, default_val=0.0, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        self.checkbox = QCheckBox(label_text)
        self.checkbox.setFixedWidth(LABEL_WIDTH)
        self.checkbox.toggled.connect(self._on_toggle)
        layout.addWidget(self.checkbox)

        self.spinbox = QDoubleSpinBox()
        self.spinbox.setSuffix(f"  {suffix}" if suffix else "")
        self.spinbox.setRange(min_val, max_val)
        self.spinbox.setDecimals(decimals)
        self.spinbox.setValue(default_val)
        self.spinbox.setEnabled(False)
        self.spinbox.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        layout.addWidget(self.spinbox, 1)

    def _on_toggle(self, checked):
        self.spinbox.setEnabled(checked)

    def is_enabled(self):
        return self.checkbox.isChecked()

    def value(self):
        return self.spinbox.value() if self.checkbox.isChecked() else None

    def set_state(self, enabled, value=None):
        self.checkbox.setChecked(bool(enabled))
        if value is not None:
            self.spinbox.setValue(value)


class CheckableGroupBox:
    """
    Mixin helper: wires a QGroupBox(checkable=True) so that unchecking it
    disables a given content widget (its actual child container).
    """

    @staticmethod
    def wire(groupbox, content_widget):
        content_widget.setEnabled(groupbox.isChecked())
        groupbox.toggled.connect(content_widget.setEnabled)


class KeyRangeTable(QWidget):
    """
    Editable table of {name: (min, max)} entries, used for the Starlight
    population/SFR/FC/BB/ret-mass age-limit dictionaries.

    Values accept plain floats or scientific notation strings (e.g. "1E7"),
    which are parsed with float() when read back.
    """

    def __init__(self, name_header="Name", min_header="Min", max_header="Max", parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels([name_header, min_header, max_header])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(90)
        self.table.setMaximumHeight(180)
        layout.addWidget(self.table)

        btn_row = QWidget()
        btn_lay = QHBoxLayout(btn_row)
        btn_lay.setContentsMargins(0, 0, 0, 0)
        btn_lay.setSpacing(6)

        btn_add = QPushButton("+ Add Row")
        btn_add.setObjectName("minorBtn")
        btn_add.setCursor(Qt.PointingHandCursor)
        btn_add.clicked.connect(self.add_row)
        btn_lay.addWidget(btn_add)

        btn_del = QPushButton("- Remove Selected")
        btn_del.setObjectName("minorBtn")
        btn_del.setCursor(Qt.PointingHandCursor)
        btn_del.clicked.connect(self._remove_selected)
        btn_lay.addWidget(btn_del)
        btn_lay.addStretch()

        layout.addWidget(btn_row)

    def add_row(self, name="", min_val="", max_val=""):
        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setItem(row, 0, QTableWidgetItem(str(name)))
        self.table.setItem(row, 1, QTableWidgetItem(str(min_val)))
        self.table.setItem(row, 2, QTableWidgetItem(str(max_val)))

    def _remove_selected(self):
        rows = sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True)
        for row in rows:
            self.table.removeRow(row)

    def set_rows(self, entries):
        """entries: iterable of (name, min, max)"""
        self.table.setRowCount(0)
        for name, min_val, max_val in entries:
            self.add_row(name, min_val, max_val)

    def get_dict(self):
        """Returns {name: (min, max)} skipping incomplete rows. Raises ValueError on bad numbers."""
        result = {}
        for row in range(self.table.rowCount()):
            name_item = self.table.item(row, 0)
            min_item = self.table.item(row, 1)
            max_item = self.table.item(row, 2)
            name = name_item.text().strip() if name_item else ""
            if not name:
                continue
            min_txt = min_item.text().strip() if min_item else ""
            max_txt = max_item.text().strip() if max_item else ""
            try:
                min_val = float(min_txt)
                max_val = float(max_txt)
            except ValueError:
                raise ValueError(
                    f"Row '{name}': min/max must be numeric (got '{min_txt}', '{max_txt}')."
                )
            result[name] = (min_val, max_val)
        return result

    def set_dict(self, d):
        entries = [(name, vals[0], vals[1]) for name, vals in (d or {}).items()]
        self.set_rows(entries)


class TargetsTableWidget(QWidget):
    """
    Editable spreadsheet-like view of the targets CSV consumed by
    Urutau.read_csv(): first column identifies each target file, any
    extra column (redshift, ebv, galaxy distance, ...) overrides that
    parameter per-target. Can be created by hand, edited, or loaded
    from/saved to an existing CSV file.
    """

    csv_path_changed = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self.table = QTableWidget(0, 1)
        self.table.setHorizontalHeaderLabels(["target"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.verticalHeader().setVisible(False)
        self.table.setMinimumHeight(140)
        self.table.setMaximumHeight(260)
        layout.addWidget(self.table)

        info = QLabel(
            "First column = target filename (relative to Targets Directory). "
            "Extra columns (e.g. redshift, ebv, galaxy distance, mask file) override that "
            "parameter per-target across every module that uses it — e.g. a 'mask file' "
            "column gives each galaxy its own STARLIGHT emission-line mask."
        )
        info.setProperty("muted", "true")
        info.setWordWrap(True)
        layout.addWidget(info)

        btn_row = QWidget()
        btn_lay = QHBoxLayout(btn_row)
        btn_lay.setContentsMargins(0, 0, 0, 0)
        btn_lay.setSpacing(6)

        def _mk(text, slot):
            btn = QPushButton(text)
            btn.setObjectName("minorBtn")
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(slot)
            btn_lay.addWidget(btn)
            return btn

        _mk("+ Row", lambda: self.add_row())
        _mk("- Row", self._remove_selected_rows)
        _mk("+ Column", self._add_column_dialog)
        _mk("- Column", self._remove_selected_column)
        btn_lay.addSpacing(12)
        _mk("Load CSV…", self._on_load_csv)
        _mk("Save CSV As…", self._on_save_csv_as)
        btn_lay.addStretch()
        layout.addWidget(btn_row)

    # -- row / column editing -------------------------------------------------

    def add_row(self, values=None):
        row = self.table.rowCount()
        self.table.insertRow(row)
        ncols = self.table.columnCount()
        values = values or []
        for c in range(ncols):
            text = str(values[c]) if c < len(values) and values[c] is not None else ""
            self.table.setItem(row, c, QTableWidgetItem(text))

    def _remove_selected_rows(self):
        rows = sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True)
        for row in rows:
            self.table.removeRow(row)

    #: Suggested per-target override keys — these match module config keys
    #: literally (Urutau matches CSV columns to config keys by exact name),
    #: e.g. "mask file" overrides StarlightOnUrutau's per-target emission
    #: line mask. Pick one or type a custom key/value pair.
    COMMON_COLUMNS = [
        "redshift", "ebv", "galaxy distance", "mask file",
        "sigma_ini", "sigma_fin", "input value", "output value",
        "timeout minutes",
    ]

    def _add_column_dialog(self):
        name, ok = QInputDialog.getItem(
            self, "Add Column",
            "Column name — pick a common override key or type your own\n"
            "(must match a module's parameter name exactly to take effect):",
            self.COMMON_COLUMNS, 0, True
        )
        name = name.strip()
        if not ok or not name:
            return
        headers = self._headers()
        if name in headers:
            QMessageBox.warning(self, "Duplicate Column", f"Column '{name}' already exists.")
            return
        col = self.table.columnCount()
        self.table.insertColumn(col)
        headers.append(name)
        self.table.setHorizontalHeaderLabels(headers)
        for row in range(self.table.rowCount()):
            self.table.setItem(row, col, QTableWidgetItem(""))

    def _remove_selected_column(self):
        cols = sorted({idx.column() for idx in self.table.selectedIndexes()}, reverse=True)
        for col in cols:
            if col == 0:
                QMessageBox.warning(self, "Cannot Remove", "The first column ('target') can't be removed.")
                continue
            self.table.removeColumn(col)

    def _headers(self):
        return [
            self.table.horizontalHeaderItem(c).text() if self.table.horizontalHeaderItem(c) else f"col{c}"
            for c in range(self.table.columnCount())
        ]

    # -- data access ------------------------------------------------------------

    def get_columns_rows(self):
        headers = self._headers()
        rows = []
        for r in range(self.table.rowCount()):
            row_vals = []
            for c in range(len(headers)):
                item = self.table.item(r, c)
                row_vals.append(item.text().strip() if item else "")
            if any(v for v in row_vals):
                rows.append(row_vals)
        return headers, rows

    def set_columns_rows(self, columns, rows):
        columns = list(columns) if columns else ["target"]
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self.table.setRowCount(0)
        for row_vals in rows or []:
            self.add_row(row_vals)

    def get_dataframe(self):
        headers, rows = self.get_columns_rows()
        return pd.DataFrame(rows, columns=headers)

    def set_dataframe(self, df):
        self.set_columns_rows(list(df.columns), df.values.tolist())

    def save_csv(self, path):
        self.get_dataframe().to_csv(path, index=False)

    # -- load / save dialogs ------------------------------------------------------

    def _on_load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Targets CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if not path:
            return
        try:
            df = pd.read_csv(path, dtype=str).fillna("")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not read CSV: {e}")
            return
        self.set_dataframe(df)
        self.csv_path_changed.emit(path)

    def _on_save_csv_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Targets CSV As", "targets.csv", "CSV Files (*.csv);;All Files (*)"
        )
        if not path:
            return
        self.save_csv(path)
        self.csv_path_changed.emit(path)


class LogConsole(QTextEdit):
    """Styled, read-only text console for log output."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setReadOnly(True)
        self.setPlaceholderText("Log output will appear here…")
        self.setMinimumHeight(120)

    def append_log(self, text):
        self.append(text)
        scrollbar = self.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def clear_log(self):
        self.clear()
