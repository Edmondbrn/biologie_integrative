from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QPlainTextEdit, QComboBox, QFileDialog, QProgressBar, QCheckBox
)
import pandas as pd
import os
import traceback

from .manual_distances_window import ManualDistancesWindow
from .app_utils import show_alert, load_stylesheet

from ..Back.distances import Distances
from ..Back.DistanceWorkerAll import DistancesWorkerAll, ParallelDistancesWorkerAll
from ..Back.distances_utils import FilterDataProt

from ..GLOBAL import *

class AllSplicingDistancesWindow(ManualDistancesWindow):
    
    def __init__(self, splice_type : str, reference_file, genomic_file):
        super().__init__(reference_file, genomic_file)
        self.splice = splice_type
        self.couple = []
        self.setWindowTitle("All Splicing Dis Calculation")

    def create_second_file_section(self):
        """
        Override this mother method to let the user selcect a folder insted of a file
        """
        group_second = QGroupBox("Folder selection")
        second_layout = QVBoxLayout(group_second)

        self.label_instruction_2 = QLabel("Please select the folder containing all the files to compare")
        self.label_instruction_2.setAlignment(Qt.AlignmentFlag.AlignCenter)
        second_layout.addWidget(self.label_instruction_2)

        file2_layout = QHBoxLayout()
        self.input_folder = QLabel("No folder selected")
        self.button_2 = QPushButton("Select input folder")
        self.button_2.clicked.connect(lambda: self.select_input_directory("Input folder"))
        file2_layout.addWidget(self.input_folder)
        file2_layout.addWidget(self.button_2)
        second_layout.addLayout(file2_layout)

        self.second_separator_label = QLabel("Please select the separator used in the files")
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

    def select_input_directory(self, dir_path):
        """
        Override this mother method to let the user select a folder instead of a file and to set the second separator label and combo box visible
        """
        file_dialog = QFileDialog()
        file_dialog.setOption(QFileDialog.Option.ShowDirsOnly)
        file_dialog.setFileMode(QFileDialog.FileMode.Directory)
        directory = file_dialog.getExistingDirectory(self, dir_path)
        if directory:
            self.input_folder.setText(f"{dir_path} : {directory}")
            self.file_dict["second"] = directory
            self.second_separator_label.setVisible(True)
            self.second_separator_combo.setVisible(True)

    def validate_files(self):
        """
        Method to validate the selected files and proceed to the pairwise selection for distance computation.
        """
        if not self.file_dict["reference"] or not self.file_dict["second"]:
            show_alert("Error", "Reference file and input directory must be selected.")
            return
        try:
            if self.file_dict["reference"].endswith(".csv"):
                self.df_ref = pd.read_csv(self.file_dict["reference"], sep=self.first_separator_combo.currentText().split(" | ")[1], engine = "python")
            if self.file_dict["reference"].endswith(".tsv"):
                self.df_ref = pd.read_csv(self.file_dict["reference"], sep="\t", engine = "python")
            if self.file_dict["reference"].endswith(".xlsx"):
                self.df_ref = pd.read_excel(self.file_dict["reference"])

            # reading the files inside the input folder
            self.dict_splicing_files = {}
            for files in os.listdir(self.file_dict["second"]):
                self.dict_splicing_files[files.split(".")[0]] = pd.read_csv(os.path.join(self.file_dict["second"], files), sep=self.second_separator_combo.currentText().split(" | ")[1], engine = "python")
                
            self.validate_button.setVisible(False) # enlève le bouton "Valider"
            self.show_column_selection()
        except Exception as e:
            show_alert("Error", f"Failed to read files: {e}")
            return

    def addColumnsSelection(self):
        """
        Method to add the column selection widgets for comparing columns from the reference dataframe and all the otehrs.
        """
        self.column_selection_label = QLabel("Select columns to compare:")
        self.group_layout.addWidget(self.column_selection_label)

        # Zone de texte pour afficher les paires
        self.comparison_label = QLabel("Comparison pairs:")
        self.comparison_text = QPlainTextEdit()
        self.comparison_text.setReadOnly(True)
        self.group_layout.addWidget(self.comparison_text)

        self.button_compare_box = QHBoxLayout()
        # Bouton "Compare"
        self.compare_button = QPushButton("Compare")
        self.compare_button.clicked.connect(self.compare_columns)
        self.button_compare_box.addWidget(self.compare_button)

        # Fixer les couples de colonnes à comparer
        try :
            for couple in self.GenerateCouple():
                pair = f"{couple[0]} - {couple[1]}"
                self.comparison_text.appendPlainText(pair)
                self.compare_pairs.append(pair) # liste pour la classe mère pour qu'elle le passe au worker
        except Exception as e:
            show_alert("Error", f"Failed to read files. Please check the column names: {e}")
            return
        self.group_layout.addLayout(self.button_compare_box)

    def GenerateCouple(self):
        """
        Method to generate the couple of columns to compare.
        Note : this method need to respect the header of the files defined in filteredRmats
        """
        self.dict_couple = {"A5SS_+": A5SS_COL, "A5SS_-": A5SS_COL, "RI_+": RI_COL, 
                            "RI_-": RI_COL, "A3SS_+": A3SS_COL, "A3SS_-": A3SS_COL, 
                            "SE_+": SE_COL, "SE_-": SE_COL, "MXE_+": MXE_COL, "MXE_-": MXE_COL}
        self.dict_splice_couples = {} # dictionnaire pour lier les couples de colonnes à chaque type de splicing
        # browse all the files in the input directory
        for splice, _ in self.dict_splicing_files.items():
            # get the columns of the reference dataframe
            for col_name_ref in REF_COUPLE:
                # get the columns of the splicing dataframe
                for col_name in self.dict_couple[f"{splice}"]:
                    self.couple.append((col_name_ref, col_name))
                    if splice not in self.dict_splice_couples:
                        self.dict_splice_couples[splice] = [(col_name_ref, col_name)]
                    else:
                        self.dict_splice_couples[splice].append((col_name_ref, col_name))
        return self.couple
    
    def addThreadsSelection(self):
        """
        Method to add the number of processing to use durinf the calculation.
        """
        box_multi = QHBoxLayout()
        self.choose_parallelisation = QLabel("Activate multithreading ? (not recommended on Windows)")
        self.choice = QCheckBox()
        self.choice.setStyleSheet(load_stylesheet(QSS_PATH))


        # fill the different boxes and add them to the layout
        box_multi.addWidget(self.choose_parallelisation)
        box_multi.addWidget(self.choice)

        self.group_layout.addLayout(box_multi)
    

    
    def compare_columns(self):
        """
        Compare toutes les paires stockées dans self.compare_pairs.
        """
        # Get all the pairs in a correct format for the Distances class
        # TODO prendre en compte l'organisme et la version de ensembl ICI
        if self.comparison_text.toPlainText() != "":
            if self.choice.isChecked():
                self.startParallelCalculation()
            else:
                self.startCalculation()
        else:
            show_alert("Error", "No pair to compare")

    def startCalculation(self):
        # Création du thread
        try:
            if self.df_ref.get("ensembl_id") is None:
                raise Exception("The 'ensembl_id' column is not found in the reference file.")
            self.worker = DistancesWorkerAll(df_ref = self.df_ref, 
                                         input_df = self.dict_splicing_files, 
                                         comparison_couples = self.dict_splice_couples,
                                         output_dir = self.output_directory.text().split(":")[1][1:], 
                                         release = RELEASE,
                                         species = SPECY,
                                         file_basename = self.file_name_space.toPlainText())
        
            self.worker.progress_changed.connect(self.updateProgressBar)
            self.worker.finished_signal.connect(self.onCalculationFinished)
            self.worker.error_signal.connect(self.onWorkerError)

            self.addProgressBar() # ajout de la barre de progression avant de lancer le calcul

            self.worker.start()
        except Exception as e:
            show_alert("Error", f"Failed in calculation\n {traceback.format_exc()}.")
            return
    
    def onWorkerError(self, error_message):
        show_alert("Error", f"Failed in calculation.\n  {error_message}.")

    def startParallelCalculation(self):
        """
        Method to initiate the parallel calculation of the distances. and to link the signals to the GUI.
        """
        try :
            if self.df_ref.get("ensembl_id") is None:
                raise Exception("The 'ensembl_id' column is not found in the reference file.")
            self.worker = ParallelDistancesWorkerAll(df_ref=self.df_ref,
                                                input_dfs=self.dict_splicing_files,
                                                comparison_couples=self.dict_splice_couples,
                                                n_processes=len(self.dict_splicing_files),
                                                release=RELEASE,
                                                species=SPECY,
                                                output_dir = self.output_directory.text().split(":")[1][1:],
                                                file_basename=self.file_name_space.toPlainText())

            self.worker.progress_changed.connect(self.updateParallelProgressBar)
            self.worker.finished_signal.connect(self.onCalculationFinished)

            self.addProgressBar()

            self.worker.start()
        except Exception as e:
            show_alert("Error", f"Failed in parallel calculation\n {traceback.format_exc()}.")
            return


    def addProgressBar(self):
        # Crée la barre de progression
        self.progress = QProgressBar()
        num_splicing_files = len(self.dict_splicing_files)
        total = len(FilterDataProt(self.df_ref)) * num_splicing_files
        self.progress.setRange(0, total)
        self.progress.setFixedWidth(300)  # Largeur fixe
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
    