import sys
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QPlainTextEdit, QDialog, QFileDialog, QComboBox, QCheckBox,
    QSpinBox, QProgressBar
)
from .app_utils import load_stylesheet, show_alert

import os
import pandas as pd
import numpy as np

from ..Back.distances import Distances
from ..Back.distances_utils import FilterDataProt
from ..Back.DistanceWorker import DistancesWorker, ParallelDistancesWorker

from ..GLOBAL import *

class ManualDistancesWindow(QDialog):
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
        self.create_reference_file_section()
        # ====================== SECTION FICHIER 2 =======================
        self.create_second_file_section()

        # ====================== SECTION SELECT OUTPUT DIRECTORY AND RESULT FILE NAMEs =======================
        self.create_output_section()

        # ====================== SECTION BOUTON VALIDATION =======================
        # Bouton de validation
        self.validate_button = QPushButton("Validate")
        self.validate_button.clicked.connect(self.validate_files)
        self.main_layout.addWidget(self.validate_button)

        self.setLayout(self.main_layout)

    def create_reference_file_section(self):
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

    def create_second_file_section(self):
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

    def create_output_section(self):
        group_output = QGroupBox("Output")
        third_layout = QVBoxLayout(group_output)
        self.label_instruction_3 = QLabel("Please select the output directory")
        self.label_instruction_3.setAlignment(Qt.AlignmentFlag.AlignCenter)
        third_layout.addWidget(self.label_instruction_3)

        output_box = QHBoxLayout()
        self.output_directory = QLabel("No directory selected")
        self.button_output = QPushButton("Select output directory")
        self.button_output.clicked.connect(lambda: self.select_output_directory("Output directory"))
        output_box.addWidget(self.output_directory)
        output_box.addWidget(self.button_output)
        third_layout.addLayout(output_box)

        self.file_name_space = QPlainTextEdit()
        self.file_name_space.setPlaceholderText("Enter the name of the result file")
        third_layout.addWidget(self.file_name_space)

        self.main_layout.addWidget(group_output)


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
                if file_path.endswith(".csv"):
                    self.first_separator_label.setVisible(True)
                    self.first_separator_combo.setVisible(True)
                self.file_dict["reference"] = file_path
            elif file_number == 2: # deuxième fichier contenant les distances
                if file_path.endswith(".csv"):
                    self.second_separator_label.setVisible(True)
                    self.second_separator_combo.setVisible(True)
                self.file_dict["second"] = file_path

    def validate_files(self):
        """
        Method to validate the selected files and proceed to the pairwise selection for distance computation.
        """
        if not self.file_dict["reference"] or not self.file_dict["second"]:
            show_alert("Error", "Both files must be selected.")
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
            show_alert("Error", f"Failed to read files: {e}")
            return
        
    def addThreadsSelection(self):
        """
        Method to add the number of processing to use durinf the calculation.
        """
        box_multi = QHBoxLayout()
        self.choose_parallelisation = QLabel("Activate multithreading ? (not recommended on Windows)")
        self.choice = QCheckBox()
        self.choice.setStyleSheet(load_stylesheet(QSS_PATH))
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
        self.group_layout.addLayout(box_multi)
        self.group_layout.addLayout(box_spin)

    def addColumnsSelection(self):
        """
        Method to add the column selection widgets for comparing columns from the two dataframes.
        """
        self.column_selection_label = QLabel("Select columns to compare:")
        self.group_layout.addWidget(self.column_selection_label)

        self.column_combo_ref = QComboBox()
        # ajout les colonnes numériques du dataframe de référence
        self.column_combo_ref.addItems(lambda: self.df_ref.select_dtypes(include=[np.number]).columns)
        self.group_layout.addWidget(self.column_combo_ref)

        self.column_combo_second = QComboBox()
        self.column_combo_second.addItems(self.df_second.select_dtypes(include=[np.number]).columns)
        self.group_layout.addWidget(self.column_combo_second)

        # Zone de texte pour afficher les paires
        self.comparison_label = QLabel("Comparison pairs:")
        self.comparison_text = QPlainTextEdit()
        self.comparison_text.setReadOnly(True)
        self.group_layout.addWidget(self.comparison_text)
        self.group_layout.addWidget(self.comparison_text)

        self.button_compare_box = QHBoxLayout()

        # Bouton d'ajout de comparaison
        self.add_comparison_button = QPushButton("Add comparison")
        self.add_comparison_button.clicked.connect(self.add_comparison)
        self.button_compare_box.addWidget(self.add_comparison_button)

        # Bouton "Compare"
        self.compare_button = QPushButton("Compare")
        self.compare_button.clicked.connect(self.compare_columns)
        self.button_compare_box.addWidget(self.compare_button)

        self.group_layout.addLayout(self.button_compare_box)

    def addProgressBar(self):
        # Crée la barre de progression
        self.progress = QProgressBar()
        self.progress.setRange(0, len(FilterDataProt(self.df_ref)))
        self.progress.setFixedWidth(300)  # Largeur fixe
        self.progress.setFormat("%p%")
        self.first_update = True
        
        # Créer un sous-layout horizontal pour centrer la barre
        self.group_progress = QGroupBox("Progress")
        self.progress_layout = QVBoxLayout(self.group_progress)
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(self.progress, alignment=Qt.AlignmentFlag.AlignCenter)
        hbox.addStretch(1)
        self.progress_layout.addLayout(hbox)

        self.layout().addWidget(self.group_progress)
        
    def show_column_selection(self):
        """
        Method to display column selection widgets for comparing columns from the two dataframes.
        """
        self.compare_pairs = []  # Liste pour stocker les paires de colonnes

        self.group_compare = QGroupBox("Column comparison")
        self.group_layout = QVBoxLayout(self.group_compare)
        # ====================== Section for the multithreading layout  =======================
        self.addThreadsSelection()
        # ====================== Section for the column selection layout  =======================
        self.addColumnsSelection()

        self.layout().addWidget(self.group_compare)


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
            show_alert("Warning", f"The pair '{pair}' already exists.")

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
            try:
                self.startParallelCalculation(comparison_list, dist.bdd)
            except Exception as e:
                show_alert("Error", "Failed to start parallel calculation.\n", e)
        else:
            try :
                self.startCalculation(comparison_list, dist.bdd)
            except Exception as e:
                show_alert("Error", "Failed to start calculation.\n", e)

        
    def startCalculation(self, comparison_list, bdd, splice_name : str = "", cpt : int = 0):
        # Création du thread
        try:
            self.worker = DistancesWorker(df_ref = self.df_ref, 
                                    df_second = self.df_second, 
                                    comparison_couples = comparison_list,
                                    output_dir = self.output_directory.text().split(":")[1][1:], 
                                    bdd = bdd,
                                    file_basename = self.file_name_space.toPlainText() + "_" +splice_name)
            
            self.worker.progress_changed.connect(self.updateProgressBar)
            self.worker.finished_signal.connect(self.onCalculationFinished)

            self.addProgressBar() # ajout de la barre de progression avant de lancer le calcul

            self.worker.start()
        except Exception as e:
            show_alert("Error", "Failed to start calculation.\n", e)
            return

    def startParallelCalculation(self, comparison_list, bdd, splice_name : str = ""):
        """
        Method to initiate the parallel calculation of the distances. and to link the signals to the GUI.
        """
        try:
            self.worker = ParallelDistancesWorker(df_ref=self.df_ref,
                                                df_splicing=self.df_second,
                                                comparison_couples=comparison_list,
                                                n_processes=self.thread_counter.value(),
                                                bdd=bdd,
                                                output_dir = self.output_directory.text().split(":")[1][1:],
                                                file_basename=self.file_name_space.toPlainText() + "_" + splice_name)

            self.worker.progress_changed.connect(self.updateParallelProgressBar)
            self.worker.finished_signal.connect(self.onCalculationFinished)

            self.addProgressBar()

            self.worker.start()
        except Exception as e:
            show_alert("Error", "Failed to start parallel calculation.\n", e)
            return

    def updateParallelProgressBar(self, rows_done: int):
        current_value = self.progress.value()
        self.progress.setValue(current_value + rows_done)

    def updateProgressBar(self, value):
        """
        Method which is solicited when an update emit is emitted by the worker.
        """
        if self.progress:
            if self.first_update:
                self.first_update = False
                self.progress.setFormat("%p%")
            self.progress.setValue(value)

    def onCalculationFinished(self):
        """
        Method which is solicited when the worker has finished its job.
        """
        # Ferme la barre de progression ou autre
        if self.progress:
            self.progress.close()
            self.layout().removeWidget(self.group_progress)
            self.group_progress.deleteLater()
            self.group_progress = None
            show_alert("Info", "Calculation finished")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = ManualDistancesWindow()
    dialog.exec()