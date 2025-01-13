from PyQt6.QtWidgets import QMainWindow, QApplication, QFileDialog, QPushButton, QVBoxLayout, QWidget, QLabel
from ..Back.parsing_rmats import getRmatsFiles, filterStrand, chooseParsing
from .app_utils import load_stylesheet, show_alert
from ..GLOBAL import QSS_PATH
import os

class ParsingRmats(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rMATS Parser")
        self.setGeometry(100, 100, 400, 300)
        self.setStyleSheet(load_stylesheet(QSS_PATH))

        layout = QVBoxLayout()

        self.rmats_label = QLabel("Select rMATS directory:")
        layout.addWidget(self.rmats_label)

        self.rmats_button = QPushButton("Choose rMATS Directory")
        self.rmats_button.clicked.connect(self.chooseRmatsDir)
        layout.addWidget(self.rmats_button)

        self.output_label = QLabel("Select Output directory:")
        layout.addWidget(self.output_label)

        self.output_button = QPushButton("Choose Output Directory")
        self.output_button.clicked.connect(self.chooseOutputDir)
        layout.addWidget(self.output_button)

        self.process_button = QPushButton("Process rMATS Files")
        self.process_button.clicked.connect(self.startProcessing)
        layout.addWidget(self.process_button)

        self.rmats_dir = ""
        self.output_dir = ""

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def chooseRmatsDir(self):
        self.rmats_dir = QFileDialog.getExistingDirectory(self, "Select rMATS Directory")
        self.rmats_label.setText(f"rMATS Directory: {self.rmats_dir}")

    def chooseOutputDir(self):
        self.output_dir = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        self.output_label.setText(f"Output Directory: {self.output_dir}")

    def startProcessing(self):
        if self.rmats_dir and self.output_dir:
            try:
                self.processRmats(self.rmats_dir, self.output_dir)
                show_alert("Info", "rMATS files processed successfully!")
            except Exception as e:
                show_alert("Error", f"An error occurred: {str(e)}")
        else:
            show_alert("Warning", "Please select both rMATS and output directories!")

    def processRmats(self, rmats_dir: str, output_dir: str):
        # Change the current working directory to the output directory
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)

        # Call the rMATS processing functions
        liste_imp = getRmatsFiles(rmats_dir)
        rmats_dico = filterStrand(liste_imp)
        chooseParsing(rmats_dico)

if __name__ == "__main__":
    app = QApplication([])
    window = ParsingRmats()
    window.show()
    app.exec()
