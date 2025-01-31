import os
import sys
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QFileDialog, QPlainTextEdit, QCheckBox, QInputDialog
)

from io import StringIO
from PyQt6.QtGui import QIcon
from .app_utils import show_alert, load_stylesheet
from ..Back.SequenceFinder import SequenceFinder
from ..Back.SequenceFinder_lazy import SequenceFinderLazy
from pandas import read_csv
from ..GLOBAL import QSS_PATH

class RNAtoDNAWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RNA to DNA Converter")
        self.setStyleSheet(load_stylesheet(QSS_PATH))
        self.resize(600, 400)
        self.setWindowIcon(QIcon("path_to_icon.png"))

        self.input_file = None
        self.output_directory = None

        # Redirection des sorties (print)
        self.stdout_buffer = StringIO()
        sys.stdout = self.stdout_buffer

        # Layout principal
        self.main_layout = QVBoxLayout(self)

        # Sections
        self.create_input_file_section()
        self.create_output_section()
        self.create_validation_button()

        # Ajout au layout principal
        self.setLayout(self.main_layout)

    def create_input_file_section(self):
        """
        Section pour sélectionner le fichier d'entrée.
        """
        group_input = QGroupBox()
        input_layout = QVBoxLayout(group_input)

        instruction_label = QLabel("Please select the RNA file to convert to DNA")
        instruction_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        input_layout.addWidget(instruction_label)

        file_layout = QHBoxLayout()
        self.input_file_label = QLabel("No file selected")
        open_file_button = QPushButton("Select File")
        open_file_button.clicked.connect(self.select_input_file)
        file_layout.addWidget(self.input_file_label)
        file_layout.addWidget(open_file_button)

        input_layout.addLayout(file_layout)

        self.choice_program = QCheckBox("Skip sequence alignement ?")
        input_layout.addWidget(self.choice_program)

        self.main_layout.addWidget(group_input)

    def select_input_file(self):
        """
        Ouvre une boîte de dialogue pour sélectionner le fichier d'entrée.
        """
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "Select RNA File", "", "CSV Files (*.csv);;TSV Files (*.tsv*);;Excel Files (*.xlsx)")
        if file_path:
            self.input_file = file_path
            self.input_file_label.setText(f"Selected: {os.path.basename(file_path)}")

    def create_output_section(self):
        """
        Section pour sélectionner le répertoire de sortie et le nom du fichier.
        """
        group_output = QGroupBox()
        output_layout = QVBoxLayout(group_output)

        instruction_label = QLabel("Please select the output directory and file name")
        instruction_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        output_layout.addWidget(instruction_label)

        directory_layout = QHBoxLayout()
        self.output_directory_label = QLabel("No directory selected")
        select_dir_button = QPushButton("Select Directory")
        select_dir_button.clicked.connect(self.select_output_directory)
        directory_layout.addWidget(self.output_directory_label)
        directory_layout.addWidget(select_dir_button)

        output_layout.addLayout(directory_layout)

        self.output_file_name = QPlainTextEdit()
        self.output_file_name.setPlaceholderText("Enter the name of the output file (e.g., converted_file.txt)")
        output_layout.addWidget(self.output_file_name)

        self.main_layout.addWidget(group_output)
    

    def select_output_directory(self):
        """
        Ouvre une boîte de dialogue pour sélectionner le répertoire de sortie.
        """
        file_dialog = QFileDialog()
        file_dialog.setOption(QFileDialog.Option.ShowDirsOnly)
        directory = file_dialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_directory = directory
            self.output_directory_label.setText(f"Selected: {directory}")

    def create_validation_button(self):
        """
        Bouton pour valider la sélection et lancer la conversion.
        """
        validate_button = QPushButton("Convert RNA to DNA")
        validate_button.clicked.connect(self.start_conversion)
        self.main_layout.addWidget(validate_button)

    def start_conversion(self):
        """
        Lancer la conversion via la méthode `start()` de `SequenceFinder`.
        """
        if not self.input_file:
            show_alert("Error", "Please select an input file.")
            return

        if not self.output_directory:
            show_alert("Error", "Please select an output directory.")
            return

        output_name = self.output_file_name.toPlainText().strip()
        if not output_name:
            show_alert("Error", "Please enter a valid output file name.")
            return

        try:
            # Lire le fichier d'entrée contenant les données ARN
            df_rna = read_csv(self.input_file, sep="\t")

            # Initialiser l'objet SequenceFinder
            if self.choice_program.isChecked():
                seq_finder = SequenceFinderLazy(data_prot=df_rna)
                seq_finder.start()  # Commence la conversion
            else:
                seq_finder = SequenceFinder(data_prot=df_rna)
                seq_finder.start()
            # Enregistrez le fichier résultant dans le répertoire de sortie
            output_path = os.path.join(self.output_directory, output_name)
            seq_finder.getDataProt().to_csv(f"{output_path}.tsv", sep="\t", index=False)

            show_alert("Info", f"RNA has been successfully converted to DNA.\nOutput saved at: {output_path}")
        except Exception as e:
            show_alert("Error", f"An error occurred during conversion: {str(e)}")

 