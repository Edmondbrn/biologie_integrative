import sys
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QPlainTextEdit, QDialog, QFileDialog, QComboBox, QCheckBox,
    QSpinBox, QProgressBar
)
from PyQt6.QtGui import QAction, QIcon
from .app_utils import load_stylesheet, show_alert

import os
import pandas as pd
import traceback
import csv

from ..Back.distances_utils import FilterDataProt
from ..Back.DistanceWorker import DistancesWorker, ParallelDistancesWorker

from ..GLOBAL import *

class ManualDistancesWindow(QDialog):
    data_signal = pyqtSignal(list)
    name_signal = pyqtSignal(list)
    def __init__(self, reference_file, genomic_file):
        super().__init__()
        self.setWindowTitle("Manual Calculation")
        self.resize(WINDOW_HEIGHT // 2, WINDOW_WIDTH // 2)
        self.df_ref, self.df_second = reference_file, genomic_file
        self.file_dict = {"reference": None, "second": None}
        self.setWindowIcon(QIcon(f"{ICON_PATH}BI_logo.png"))
        self.species, self.release = self.release_reader()
        self.release = int(self.release)

        self.main_layout = QVBoxLayout(self)

        # Titre
        self.label_title = QLabel("Welcome in manual distance calculation home")
        self.label_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label_title.setStyleSheet("font-size: 20px; font-weight: bold;")
        self.main_layout.addWidget(self.label_title)

        # ====================== SECTION FICHIER 1 =======================
        if reference_file is not None:
            pass
        else:
            self.create_reference_file_section()
        # ====================== SECTION FICHIER 2 =======================
        if genomic_file is not None:
            pass
        else:
            self.create_second_file_section()

        # ====================== SECTION SELECT OUTPUT DIRECTORY AND RESULT FILE NAMEs =======================
        self.create_output_section()

        # ====================== SECTION BOUTON VALIDATION =======================
        # Bouton de validation
        self.validate_button = QPushButton("Validate")
        self.validate_button.clicked.connect(self.validate_files)
        self.validate_button.clicked.connect(lambda: self.send_files(self.df_ref, self.df_second))

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
                self.file_dict["reference"] = file_path
            elif file_number == 2: # deuxième fichier contenant les distances
                self.file_dict["second"] = file_path

    def validate_files(self):
        """
        Method to validate the selected files and proceed to the pairwise selection for distance computation.
        """
        if self.df_ref is None or self.df_second is None:
            if not self.file_dict["reference"] and not self.file_dict["second"]:
                show_alert("Error", "Both files must be selected.")
                return
            try:
                if self.file_dict["reference"].endswith(".csv"):
                    sep = self.detect_separator(self.file_dict["reference"])
                    index = self.index_column_detector(self.file_dict["reference"], sep)
                    self.df_ref = pd.read_csv(self.file_dict["reference"], sep=sep, index_col=index, engine = "python")
                if self.file_dict["second"].endswith(".csv"):
                    sep = self.detect_separator(self.file_dict["second"])
                    index = self.index_column_detector(self.file_dict["second"], sep)
                    self.df_second = pd.read_csv(self.file_dict["second"], sep=sep, index_col=index,engine = "python")
                if self.file_dict["reference"].endswith(".tsv"):
                    self.df_ref = pd.read_csv(self.file_dict["reference"], sep="\t", engine = "python")
                if self.file_dict["second"].endswith(".tsv"):
                    self.df_second = pd.read_csv(self.file_dict["second"], sep="\t", engine = "python")
                if self.file_dict["reference"].endswith(".xlsx"):
                    index = self.index_column_detector_excel(self.file_dict["reference"])
                    self.df_ref = pd.read_excel(self.file_dict["reference"], index_col=index)
                if self.file_dict["second"].endswith(".xlsx"):
                    index = self.index_column_detector_excel(self.file_dict["second"])
                    self.df_second = pd.read_excel(self.file_dict["second"], index_col=index)
                self.validate_button.setVisible(False)
                self.show_column_selection()
            except Exception as e:
                show_alert("Error", f"Failed to read files: {e}")
                return
        else:
            self.show_column_selection()
        
    def send_files(self, first_file, second_file):
        files = [first_file, second_file]
        names = [os.path.basename(self.file_dict["reference"]), os.path.basename(self.file_dict["second"])]
        self.data_signal.emit(files)
        self.name_signal.emit(names)
        
    def addThreadsSelection(self):
        """
        Method to add the number of processing to use during the calculation.
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
        self.column_combo_ref.addItems(self.df_ref.columns)
        self.group_layout.addWidget(self.column_combo_ref)

        self.column_combo_second = QComboBox()
        self.column_combo_second.addItems(self.df_second.columns)
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

        # bouton remove le dernier couple de colonnes à avoir été ajouté
        self.remove_comparison_button = QPushButton("Remove comparisons")
        self.remove_comparison_button.clicked.connect(self.clear_comparison)
        self.button_compare_box.addWidget(self.remove_comparison_button)

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

    def clear_comparison(self):
        """
        Efface la dernière paire de colonnes ajoutée.
        """
        if len(self.compare_pairs) > 0:
            self.compare_pairs.pop()
            self.comparison_text.clear()
            for pair in self.compare_pairs:
                self.comparison_text.appendPlainText(pair)
        else:
            show_alert("Warning", "No comparison pair to remove.")

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

        if self.comparison_text.toPlainText() != "":         
            if self.choice.isChecked():
                self.startParallelCalculation(comparison_list)
            else:
                self.startCalculation(comparison_list)
        else:
            show_alert("Error", "No comparison pair selected.")
            return


        
    def startCalculation(self, comparison_list, splice_name : str = ""):
        # Création du thread
        try:
            # TODO quand on doit télécharger / charger un génome, on doit chnager le svaleurs dans GLOBAL.py
            if self.df_ref.get("ensembl_id") is None:
                raise Exception("The 'ensembl_id' column is not found in the reference file.")
            self.worker = DistancesWorker(df_ref = self.df_ref, 
                                    df_second = self.df_second, 
                                    comparison_couples = comparison_list,
                                    output_dir = self.output_directory.text().split(":", maxsplit=1)[1].strip(), 
                                    release = self.release,
                                    species = self.species,
                                    file_basename = self.file_name_space.toPlainText() + "_" +splice_name)
            
            self.worker.progress_changed.connect(self.updateProgressBar)
            self.worker.finished_signal.connect(self.onCalculationFinished)
            self.worker.error_signal.connect(self.onWorkerError)
            self.addProgressBar() # ajout de la barre de progression avant de lancer le calcul

            self.worker.start()
        except Exception as e:
            show_alert("Error", f"Failed in calculation.\n  {traceback.format_exc()}.")
            return
        

    def onWorkerError(self, error_message):
        show_alert("Error", f"Failed in calculation.\n  {error_message}.")

    def startParallelCalculation(self, comparison_list, splice_name : str = ""):
        """
        Method to initiate the parallel calculation of the distances. and to link the signals to the GUI.
        """
        try:
            if self.df_ref.get("ensembl_id") is None:
                raise Exception("The 'ensembl_id' column is not found in the reference file.")
            # TODO quand on doit télécharger / charger un génome, on doit chnager le svaleurs dans GLOBAL.py
            self.worker = ParallelDistancesWorker(df_ref=self.df_ref,
                                                df_splicing=self.df_second,
                                                comparison_couples=comparison_list,
                                                n_processes=self.thread_counter.value(),
                                                release=self.release,
                                                species=self.species,
                                                output_dir = self.output_directory.text().split(":")[1][1:],
                                                file_basename=self.file_name_space.toPlainText() + "_" + splice_name)

            self.worker.progress_changed.connect(self.updateParallelProgressBar)
            self.worker.finished_signal.connect(self.onCalculationFinished)

            self.addProgressBar()

            self.worker.start()
        except Exception as e:
            show_alert("Error", f"Failed in parallel calculcation.\n {traceback.format_exc()}..")
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

    def release_reader(self, file_path):
        lines = []
        with open(file_path, "r") as file:
            for line in file:
                lines.append(line.strip())
            return lines
    
    def detect_separator(self, file_path):
        """
        Detect the separator used in a CSV file.
        """
        with open(file_path, 'r') as file:
            # Read the first line of the file
            sample = file.read(1024)
            # Try different separators
            sniffer = csv.Sniffer()
            sep = sniffer.sniff(sample).delimiter
            return sep
        return ','
    
    def index_column_detector(self, file_path, sep):
        file = pd.read_csv(file_path, sep=sep, nrows=10)
        if file.columns[0] == "" or 'Unnamed: 0' in file.columns[0]:
            return 0
        else:
            return None
        

    def index_column_detector_excel(self, file_path):
        file = pd.read_excel(file_path, nrows=10)
        if file.columns[0] == "" or 'Unnamed: 0' in file.columns[0]:
            return 0
        else:
            return None

if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = ManualDistancesWindow()
    dialog.exec()