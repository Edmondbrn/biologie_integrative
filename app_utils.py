import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QToolBar, QMessageBox, QDialog, QFileDialog, QComboBox
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt
from GLOBAL import *
import os
import pandas as pd

class FileDialogManual(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Manual Calculation")
        self.resize(WINDOW_HEIGHT // 2, WINDOW_WIDTH // 2)

        self.file_dict = {"reference": None, "second": None}

        main_layout = QVBoxLayout(self)

        # Titre
        self.label_title = QLabel("Welcome in manual distance calculation home")
        self.label_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label_title.setStyleSheet("font-size: 20px; font-weight: bold;")
        main_layout.addWidget(self.label_title)

        # ====================== SECTION FICHIER 1 =======================

        group_ref = QGroupBox("Reference file")
        ref_layout = QVBoxLayout(group_ref)

        # Instruction
        self.label_instruction_1 = QLabel("Please select the reference file to proceed")
        self.label_instruction_1.setAlignment(Qt.AlignmentFlag.AlignCenter)
        ref_layout.addWidget(self.label_instruction_1)

        # Bouton + label fichier sélectionné
        file1_layout = QHBoxLayout()
        self.first_file_path = QLabel("No file selected")
        self.button_1 = QPushButton("Open file editor")
        self.button_1.clicked.connect(lambda: self.open_file_dialog(self.first_file_path, 1))
        file1_layout.addWidget(self.first_file_path)
        file1_layout.addWidget(self.button_1)
        ref_layout.addLayout(file1_layout)

        # Label + combo pour séparateur
        self.first_separator_label = QLabel("Please select the separator used in the file")
        self.first_separator_label.setVisible(False)
        self.first_separator_combo = QComboBox()
        self.first_separator_combo.addItem("Comma | ,")
        self.first_separator_combo.addItem("Semicolon | ;")
        self.first_separator_combo.addItem("Tabulation | \\t")
        self.first_separator_combo.addItem("Space |  ")
        self.first_separator_combo.setVisible(False)
        ref_layout.addWidget(self.first_separator_label)
        ref_layout.addWidget(self.first_separator_combo)

        main_layout.addWidget(group_ref)

        # ====================== SECTION FICHIER 2 =======================
        group_second = QGroupBox("Second file")
        second_layout = QVBoxLayout(group_second)

        self.label_instruction_2 = QLabel("Please select the second file to proceed")
        self.label_instruction_2.setAlignment(Qt.AlignmentFlag.AlignCenter)
        second_layout.addWidget(self.label_instruction_2)

        file2_layout = QHBoxLayout()
        self.second_file_path = QLabel("No file selected")
        self.button_2 = QPushButton("Open file editor")
        self.button_2.clicked.connect(lambda: self.open_file_dialog(self.second_file_path, 2))
        file2_layout.addWidget(self.second_file_path)
        file2_layout.addWidget(self.button_2)
        second_layout.addLayout(file2_layout)

        self.second_separator_label = QLabel("Please select the separator used in the file")
        self.second_separator_label.setVisible(False)
        self.second_separator_combo = QComboBox()
        self.second_separator_combo.addItem("Comma | ,")
        self.second_separator_combo.addItem("Semicolon | ;")
        self.second_separator_combo.addItem("Tabulation | \\t")
        self.second_separator_combo.addItem("Space |  ")
        self.second_separator_combo.setVisible(False)
        second_layout.addWidget(self.second_separator_label)
        second_layout.addWidget(self.second_separator_combo)

        main_layout.addWidget(group_second)

        # Bouton de validation
        self.validate_button = QPushButton("Validate")
        self.validate_button.clicked.connect(self.validate_files)
        main_layout.addWidget(self.validate_button)

        self.setLayout(main_layout)

    def open_file_dialog(self, label: QLabel, file_number: int):
        """
        Method that opens a file dialog to select a file.
        There are filters for CSV, TSV and Excel files.
        """
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "Open File", "", "CSV Table (*.csv);;TSV table (*.tsv);;Excel Table (*.xlsx)")
        if file_path:
            label.setText(f"Selected file : {os.path.basename(file_path)}")
            if file_path.endswith(".csv"):
                if file_number == 1: # récupère le premier fichier de référence
                    self.first_separator_label.setVisible(True)
                    self.first_separator_combo.setVisible(True)
                    self.file_dict["reference"] = file_path
                elif file_number == 2: # deucième fichier contenant les distances
                    self.second_separator_label.setVisible(True)
                    self.second_separator_combo.setVisible(True)
                    self.file_dict["second"] = file_path

    def validate_files(self):
        """
        Method to validate the selected files and proceed to the pairwise selection for distance computation.
        """
        if not self.file_dict["reference"] or not self.file_dict["second"]:
            self.show_alert("Error", "Both files must be selected.")
            return
        try:
            df_ref = pd.read_csv(self.file_dict["reference"], sep=self.first_separator_combo.currentText().split(" | ")[1], engine = "python")
            df_second = pd.read_csv(self.file_dict["second"], sep=self.second_separator_combo.currentText().split(" | ")[1], engine = "python")
        except Exception as e:
            self.show_alert("Error", f"Failed to read files: {e}")
            return
        
    def show_alert(self, title: str, message: str):
        """
        Method to display an alert dialog with a title and a specific message.
        """
        alert = QMessageBox()
        alert.setWindowTitle(title)
        alert.setText(message)
        alert.setIcon(QMessageBox.Icon.Warning)
        alert.setStandardButtons(QMessageBox.StandardButton.Ok)
        alert.exec()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = FileDialogManual()
    dialog.exec()