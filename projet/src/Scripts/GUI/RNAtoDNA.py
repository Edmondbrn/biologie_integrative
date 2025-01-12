import os
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton,
    QLabel, QFileDialog, QPlainTextEdit, QProgressBar
)
from PyQt6.QtGui import QIcon
from .app_utils import show_alert
from ..Back.SequenceFinder import SequenceFinder
from pandas import read_csv
from ..Back.distances_utils import FilterDataProt


class RNAtoDNAWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RNA to DNA Converter")
        self.resize(600, 400)
        self.setWindowIcon(QIcon("path_to_icon.png"))

        self.input_file = None
        self.output_directory = None

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
        group_input = QGroupBox("Input File")
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
        self.main_layout.addWidget(group_input)

    def select_input_file(self):
        """
        Ouvre une boîte de dialogue pour sélectionner le fichier d'entrée.
        """
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "Select RNA File", "", "Text Files (*.txt);;All Files (*)")
        if file_path:
            self.input_file = file_path
            self.input_file_label.setText(f"Selected: {os.path.basename(file_path)}")

    def create_output_section(self):
        """
        Section pour sélectionner le répertoire de sortie et le nom du fichier.
        """
        group_output = QGroupBox("Output Directory")
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
        validate_button.clicked.connect(self.convert_rna_to_dna)
        self.main_layout.addWidget(validate_button)

    def convert_rna_to_dna(self):
        """
        Effectue la conversion du fichier RNA en fichier DNA.
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

        # Ajouter la barre de progression
        self.addProgressBar()

        # Créer un thread pour effectuer la conversion en arrière-plan
        self.worker_thread = ConversionWorker(self.input_file, self.output_directory, output_name)
        self.worker_thread.progress_signal.connect(self.updateProgressBar)
        self.worker_thread.finished.connect(self.on_conversion_finished)
        self.worker_thread.start()

    def addProgressBar(self):
        # Crée la barre de progression
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)  # Valeur maximale fixée à 100 pour commencer
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

    def updateProgressBar(self, value):
        """
        Mise à jour de la barre de progression.
        """
        if self.progress:
            if self.first_update:
                self.first_update = False
                self.progress.setFormat("%p%")
            self.progress.setValue(value)

    def on_conversion_finished(self):
        """
        Action à réaliser après la conversion (ex: informer l'utilisateur).
        """
        show_alert("Info", "RNA to DNA conversion finished!")
        self.progress.setValue(100)  # Assurez-vous que la barre de progression est remplie à 100%
        self.worker_thread.quit()

class ConversionWorker(QThread):
    progress_signal = pyqtSignal(int)

    def __init__(self, input_file, output_directory, output_name):
        super().__init__()
        self.input_file = input_file
        self.output_directory = output_directory
        self.output_name = output_name

    def run(self):
        """
        Effectue la conversion en arrière-plan.
        """
        try:
            df_rna = read_csv(self.input_file, sep="\t")
            seq_finder = SequenceFinder(data_prot=df_rna)

            # Simuler la progression (ici on suppose que la conversion se fait en étapes)
            total_steps = len(df_rna)
            for step, row in enumerate(df_rna.itertuples(), start=1):
                # Conversion (par exemple, remplace "U" par "T")
                seq_finder.convert_rna_to_dna(row)
                
                # Signal de progression
                progress = int((step / total_steps) * 100)
                self.progress_signal.emit(progress)

            # Sauvegarde du fichier converti
            output_path = os.path.join(self.output_directory, self.output_name)
            seq_finder.__data_prot.to_csv(output_path, sep="\t", index=False)

        except Exception as e:
            print(f"An error occurred: {str(e)}")
            self.progress_signal.emit(100)  # Remet la progression à 100% en cas d'erreur
