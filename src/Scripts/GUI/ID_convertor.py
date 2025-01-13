import mygene
import pandas as pd
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QProgressBar, QFileDialog, QMessageBox, QApplication, QPushButton
)
from PyQt6.QtGui import QIcon
from PyQt6.QtCore import QThread, pyqtSignal
import sys

from ..GLOBAL import *
from ..Back.Id_convertor import add_ensembl_ids

class IDConversionDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("NCBI to Ensembl Converter")
        self.setWindowIcon(QIcon(f"{ICON_PATH}address-book-blue.png"))
        self.resize(600, 200)

        self.file_path = None
        self.thread = None

        # Layout principal
        layout = QVBoxLayout(self)

        # Barre de progression
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        # Bouton de conversion
        button_convert_ID = QPushButton("Select File and Convert")
        button_convert_ID.clicked.connect(self.file_loader_ID)
        layout.addWidget(button_convert_ID)

    def file_loader_ID(self):
        """Charger un fichier et démarrer la conversion."""
        file_path, _ = QFileDialog.getOpenFileName(self, "File Explorer", "", "TSV Files (*.tsv);;CSV Files (*.csv);;All Files (*)")
        if file_path:
            self.file_path = file_path
            self.start_conversion()

    def start_conversion(self):
        """Démarrer la conversion dans un thread séparé."""
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)

        self.thread = IDConversionThread(self.file_path)
        self.thread.progress.connect(self.progress_bar.setValue)
        self.thread.finished.connect(self.on_conversion_finished)
        self.thread.error.connect(self.on_conversion_error)
        self.thread.start()

    def on_conversion_finished(self, output_path):
        """Appelé lorsque la conversion est terminée."""
        self.progress_bar.setVisible(False)
        QMessageBox.information(self, "Conversion Complete", f"Conversion finished! File saved at:\n{output_path}")

    def on_conversion_error(self, error):
        """Appelé en cas d'erreur pendant la conversion."""
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Conversion Error", f"An error occurred:\n{error}")

class IDConversionThread(QThread):
    progress = pyqtSignal(int)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path

    def run(self):
        try:
            add_ensembl_ids(self.file_path, self.progress)
            output_path = self.file_path + "_converted"
            self.finished.emit(output_path)
        except Exception as e:
            self.error.emit(str(e))

def convert_refseq_to_ensembl(refseq_ids):
    """
    Convertit une liste d'identifiants RefSeq en identifiants Ensembl
    en utilisant la librairie MyGene.
    
    refseq_ids: liste de string (e.g. ["NM_001301415", "NR_002847", ...])
    Retour: liste de dictionnaires contenant les résultats.
    """
    mg = mygene.MyGeneInfo()
    results = mg.querymany(refseq_ids, 
                           scopes='refseq', 
                           fields='ensembl.transcript', 
                           species='mouse')
    return results

def __init__(self, file_path):
        super().__init__()
        self.file_path = file_path

def run(self):
    try:
        add_ensembl_ids(self.file_path, self.progress)
        output_path = self.file_path + "_converted"
        self.finished.emit(output_path)
    except Exception as e:
        self.error.emit(str(e))

# ===================== Point d'entrée =====================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = IDConversionDialog()
    dialog.exec()
    sys.exit(app.exec())