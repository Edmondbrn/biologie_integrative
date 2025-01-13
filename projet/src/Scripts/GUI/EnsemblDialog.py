from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, 
QApplication, QMessageBox, QProgressBar)
from PyQt6.QtCore import pyqtSignal, QThread
import pyensembl as pb
from ..GLOBAL import *

class EnsemblDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Sélection de l'espèce et release Ensembl")

        layout = QVBoxLayout(self)

        # Champ pour l'espèce
        species_layout = QHBoxLayout()
        species_label = QLabel("Species :")
        self.species_input = QLineEdit()
        species_layout.addWidget(species_label)
        species_layout.addWidget(self.species_input)
        layout.addLayout(species_layout)

        # Champ pour la release
        release_layout = QHBoxLayout()
        release_label = QLabel("Release :")
        self.release_input = QLineEdit()
        release_layout.addWidget(release_label)
        release_layout.addWidget(self.release_input)
        layout.addLayout(release_layout)

        # Bouton de validation
        btn_ok = QPushButton("Validate")
        btn_ok.clicked.connect(self.on_validate)
        layout.addWidget(btn_ok)

    def on_validate(self):
        species = self.species_input.text().strip()
        release = self.release_input.text().strip()

        if not species or not release:
            QMessageBox.warning(self, "Input Error", "Both species and release fields must be filled.")
            return

        try:
            release = int(release)
        except ValueError:
            QMessageBox.warning(self, "Input Error", "Release must be an integer.")
            return

        if not self.is_valid_species(species):
            QMessageBox.warning(self, "Input Error", f"Species '{species}' is not supported.")
            return

        # Mettre à jour les variables globales
        global SPECIES, RELEASE
        SPECIES = species
        RELEASE = release

        # Créer et afficher la barre de progression
        progress_bar = QProgressBar(self)
        progress_bar.setRange(0, 100)
        progress_bar.show()

        # Lancer le téléchargement
        download_ensembl_release(species, release, progress_bar)

        self.accept()  # Ferme la boîte de dialogue

    def is_valid_species(self, species):
        # Vérifie si l'espèce est supportée par pyensembl
        try:
            pb.EnsemblRelease(release=75, species=species)  # Utilise une release existante pour vérifier
            return True
        except ValueError:
            return False

    def get_values(self):
        return self.species_input.text(), self.release_input.text()
    
    
class DownloadThread(QThread):
    progress = pyqtSignal(int)
    finished = pyqtSignal()
    error = pyqtSignal(str)

    def __init__(self, species, release):
        super().__init__()
        self.species = species
        self.release = release

    def run(self):
        try:
            ensembl_release = pb.EnsemblRelease(release=int(self.release), species=self.species)
            ensembl_release.download(progress_callback=self.update_progress)
            ensembl_release.index()
            self.finished.emit()
        except Exception as e:
            self.error.emit(str(e))

    def update_progress(self, progress):
        self.progress.emit(progress)

def ask_ensembl_info():
    dialog = EnsemblDialog()
    if dialog.exec() == QDialog.DialogCode.Accepted:
        info = dialog.get_values()
        return info
    return None, None


def download_ensembl_release(species, release, progress_bar):
    thread = DownloadThread(species, release)
    thread.progress.connect(progress_bar.setValue)
    thread.finished.connect(lambda: QMessageBox.information(None, "Download Complete", "Download and indexing complete."))
    thread.error.connect(lambda error: QMessageBox.critical(None, "Download Error", f"Failed to download Ensembl release: {error}"))
    thread.start()


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    species, release = ask_ensembl_info()
    if species and release:
        print(f"Espèce : {species}, Release : {release}")
    sys.exit(app.exec())