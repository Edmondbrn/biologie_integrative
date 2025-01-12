import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QPlainTextEdit, QMessageBox, QDialog, QFileDialog, QComboBox
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt
from GLOBAL import *
import os
import pandas as pd
from distances import Distances

class FileDialogManual(QDialog):
    def __init__(self, _file, genomic_file):
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
            if file_number == 1: # récupère le premier fichier de référence
                self.first_separator_label.setVisible(True)
                self.first_separator_combo.setVisible(True)
                self.file_dict["reference"] = file_path
            elif file_number == 2: # deuxième fichier contenant les distances
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
            if self.file_dict["reference"].endswith(".csv"):
                self.df_ref = pd.read_csv(self.file_dict["reference"], sep=self.first_separator_combo.currentText().split(" | ")[1], engine = "python")
            if self.file_dict["second"].endswith(".csv"):
                self.df_second = pd.read_csv(self.file_dict["second"], sep=self.second_separator_combo.currentText().split(" | ")[1], engine = "python")
            if self.file_dict["reference"].endswith(".tsv"):
                self.df_ref = pd.read_csv(self.file_dict["reference"], sep="\t", engine = "python")
            if self.file_dict["second"].endswith(".tsv"):
                self.df_second = pd.read_csv(self.file_dict["second"], sep="\t", engine = "python")
            if self.file_dict["reference"].endswith(".xlsx"):
                self.df_ref = pd.read_excel(self.file_dict["reference"])
            if self.file_dict["second"].endswith(".xlsx"):
                self.df_second = pd.read_excel(self.file_dict["second"])
            self.validate_button.setVisible(False)
            self.show_column_selection()
        except Exception as e:
            self.show_alert("Error", f"Failed to read files: {e}")
            return
        
    def show_column_selection(self):
        """
        Method to display column selection widgets for comparing columns from the two dataframes.
        """
        self.compare_pairs = []  # Liste pour stocker les paires de colonnes

        group_compare = QGroupBox("Column comparison")
        group_layout = QVBoxLayout(group_compare)

        self.column_selection_label = QLabel("Select columns to compare:")
        group_layout.addWidget(self.column_selection_label)

        self.column_combo_ref = QComboBox()
        self.column_combo_ref.addItems(self.df_ref.columns)
        group_layout.addWidget(self.column_combo_ref)

        self.column_combo_second = QComboBox()
        self.column_combo_second.addItems(self.df_second.columns)
        group_layout.addWidget(self.column_combo_second)

        # Zone de texte pour afficher les paires
        self.comparison_label = QLabel("Comparison pairs:")
        self.comparison_text = QPlainTextEdit()
        self.comparison_text.setReadOnly(True)
        group_layout.addWidget(self.comparison_text)
        group_layout.addWidget(self.comparison_text)

        self.button_compare_box = QHBoxLayout()

        # Bouton d'ajout de comparaison
        self.add_comparison_button = QPushButton("Add comparison")
        self.add_comparison_button.clicked.connect(self.add_comparison)
        self.button_compare_box.addWidget(self.add_comparison_button)

        # Bouton "Compare"
        self.compare_button = QPushButton("Compare")
        self.compare_button.clicked.connect(self.compare_columns)
        self.button_compare_box.addWidget(self.compare_button)

        group_layout.addLayout(self.button_compare_box)
        self.layout().addWidget(group_compare)


    def add_comparison(self):
        """
        Ajoute la paire de colonnes choisie dans une zone de texte si elle n'existe pas déjà.
        """
        col_ref = self.column_combo_ref.currentText()
        col_second = self.column_combo_second.currentText()
        pair = f"{col_ref} - {col_second}"

        if pair not in self.compare_pairs:
            self.compare_pairs.append(pair)
            self.comparison_text.appendPlainText(pair)
        else:
            self.show_alert("Info", f"The pair '{pair}' already exists.")

    def compare_columns(self):
        """
        Compare toutes les paires stockées dans self.compare_pairs.
        """
        # # Exemples d'itération des paires
        comparison_list = []
        for pair in self.compare_pairs:
            couple = tuple(pair.split(" - "))
            comparison_list.append(couple)
        dist = Distances()
        dist.start_manual(self.df_ref, self.df_second, comparison_list)

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
