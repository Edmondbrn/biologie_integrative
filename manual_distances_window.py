import sys
from PyQt6.QtWidgets import (
    QApplication, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QPlainTextEdit, QMessageBox, QDialog, QFileDialog, QComboBox, QCheckBox,
    QSpinBox, QProgressBar, QSizePolicy
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt
from GLOBAL import *
import os
import pandas as pd
from distances import Distances
from DistanceWorker import DistancesWorker
from app_utils import load_stylesheet

class FileDialogManual(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Manual Calculation")
        self.resize(WINDOW_HEIGHT // 2, WINDOW_WIDTH // 2)

        self.file_dict = {"reference": None, "second": None}

        self.main_layout = QVBoxLayout(self)

        # Titre
        self.label_title = QLabel("Welcome in manual distance calculation home")
        self.label_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label_title.setStyleSheet("font-size: 20px; font-weight: bold;")
        self.main_layout.addWidget(self.label_title)

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

        self.main_layout.addWidget(group_ref)

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

        self.main_layout.addWidget(group_second)

        # ====================== SECTION SELECT OUTPUT DIRECTORY AND RESULT FILE NAMEs =======================
        group_output = QGroupBox("Output")
        third_layout = QVBoxLayout(group_output)
        self.label_instruction_3 = QLabel("Please select the output directory")
        self.label_instruction_3.setAlignment(Qt.AlignmentFlag.AlignCenter)
        third_layout.addWidget(self.label_instruction_3)

        output_box = QHBoxLayout()
        self.output_directory = QLabel("No directory selected")
        self.button_output = QPushButton("Select output directory")
        self.button_output.clicked.connect(lambda :self.select_output_directory("Output directory"))
        output_box.addWidget(self.output_directory)
        output_box.addWidget(self.button_output)
        third_layout.addLayout(output_box)

        self.file_name_space = QPlainTextEdit()
        self.file_name_space.setPlaceholderText("Enter the name of the result file")
        third_layout.addWidget(self.file_name_space)

        self.main_layout.addWidget(group_output)

        # ====================== SECTION BOUTON VALIDATION =======================
        # Bouton de validation
        self.validate_button = QPushButton("Validate")
        self.validate_button.clicked.connect(self.validate_files)
        self.main_layout.addWidget(self.validate_button)

        self.setLayout(self.main_layout)

    def select_output_directory(self, dir_path : str):
        """
        Method to select the output directory for the results.
        """
        file_dialog = QFileDialog()
        file_dialog.setOption(QFileDialog.Option.ShowDirsOnly)
        file_dialog.setFileMode(QFileDialog.FileMode.Directory)
        directory = file_dialog.getExistingDirectory(self, dir_path)
        if directory:
            self.output_directory.setText(f"{dir_path} : {directory}")

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
        # ====================== Section for the multithreading layout  =======================
        box_multi = QHBoxLayout()
        self.choose_parallelisation = QLabel("Activate multithreading ? (not recommended on Windows)")
        self.choice = QCheckBox()
        self.choice.setStyleSheet(load_stylesheet("styles.css"))
        box_spin = QHBoxLayout()
        self.thread_counter = QSpinBox()
        self.thread_counter.setRange(2, 32)  # Définir la plage de valeurs pour le compteur de threads
        self.thread_counter.setValue(4)  # Valeur par défaut
        self.thread_counter.setVisible(False)  # Masquer le compteur de threads par défaut
        self.thread_counter.setPrefix("Threads: ")  # Préfixe pour le compteur de threads
        self.thread_counter.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Aligner le texte au centre
        self.thread_counter.setFixedWidth(100)  # Définir une largeur fixe pour le compteur de threads
        # Connecte le signal de la case à cocher à une méthode pour afficher/masquer le compteur de threads
        self.choice.stateChanged.connect(lambda: self.thread_counter.setVisible(not self.thread_counter.isVisible()))
        # fill the different boxes and add them to the layout
        box_multi.addWidget(self.choose_parallelisation)
        box_multi.addWidget(self.choice)
        box_spin.addStretch(1)
        box_spin.addWidget(self.thread_counter)
        box_spin.addStretch(1)
        group_layout.addLayout(box_multi)
        group_layout.addLayout(box_spin)
        # ====================== Section for the column selection layout  =======================
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
            self.show_alert("Warning", f"The pair '{pair}' already exists.")

    def compare_columns(self):
        """
        Compare toutes les paires stockées dans self.compare_pairs.
        """
        # Get all the pairs in a correct format for the Distances class
        comparison_list = []
        for pair in self.compare_pairs:
            couple = tuple(pair.split(" - "))
            comparison_list.append(couple)
        # TODO prendre en compte l'organisme et la version de ensembl ICI
        dist = Distances()
        if self.choice.isChecked():
            dist.parallel_start_manual(df_ref = self.df_ref, df_splicing = self.df_second, comparison_couples = comparison_list, 
                                       output_dir = self.output_directory.text().split(":")[1][1:], output_basename = self.file_name_space.toPlainText()  , 
                                       n_cores = self.thread_counter.value())
        else:
            self.startCalculation(comparison_list, dist.bdd)

        
    def startCalculation(self, comparison_list, bdd):
        # Création du thread
        self.worker = DistancesWorker(df_ref = self.df_ref, 
                                   df_second = self.df_second, 
                                   comparison_couples = comparison_list,
                                   output_dir = self.output_directory.text().split(":")[1][1:], 
                                   bdd = bdd,
                                   file_basename = self.file_name_space.toPlainText())
        
        self.worker.progress_changed.connect(self.updateProgressBar)
        self.worker.finished_signal.connect(self.onCalculationFinished)

        # Crée la barre de progression 
        # Crée la barre de progression
        self.progress = QProgressBar()
        self.progress.setRange(0, len(self.df_ref))
        self.progress.setFixedWidth(300)  # Largeur fixe

        # Ajouter la barre de progression au layout pour la centrer
        # Supposons que vous ayez un layout principal self.layout() déjà existant
        
        # Créer un sous-layout horizontal pour centrer la barre
        self.group_progress = QGroupBox("Progress")
        self.progress_layout = QVBoxLayout(self.group_progress)
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.progress, alignment=Qt.AlignmentFlag.AlignCenter)
        hbox.addStretch(1)
        self.progress_layout.addLayout(hbox)

        self.layout().addWidget(self.group_progress)
        self.worker.start()

    def updateProgressBar(self, value):
        if self.progress:
            self.progress.setValue(value)

    def onCalculationFinished(self):
        # Ferme la barre de progression ou autre
        if self.progress:
            self.progress.close()
            self.layout().removeWidget(self.group_progress)
            self.group_progress.deleteLater()
            self.group_progress = None
            self.show_alert("Info", "Calculation finished")

    def show_alert(self, title: str, message: str):
        """
        Method to display an alert dialog with a title and a specific message.
        """
        dict_logo = {"Info": QMessageBox.Icon.Information, "Error": QMessageBox.Icon.Critical, "Warning": QMessageBox.Icon.Warning}
        alert = QMessageBox()
        alert.setWindowTitle(title)
        alert.setText(message)
        alert.setIcon(dict_logo[title])
        alert.setStandardButtons(QMessageBox.StandardButton.Ok)
        alert.exec()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = FileDialogManual()
    dialog.exec()