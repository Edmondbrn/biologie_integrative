from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QApplication, QMessageBox, QProgressBar
)
from PyQt6.QtCore import pyqtSignal, QThread
from PyQt6.QtGui import QIcon
import pyensembl as pb
import sys
import os

from ..GLOBAL import *


class EnsemblDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon(f"{ICON_PATH}BI_logo.png"))
        self.setWindowTitle("Sélection de l'espèce et release Ensembl")
        self.resize(400, 200)

        # Layout principal
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

        # Barre de progression
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        # Bouton de validation
        btn_ok = QPushButton("Validate")
        btn_ok.clicked.connect(self.on_validate)
        layout.addWidget(btn_ok)

        self.download_thread = None

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
        # Écriture des informations dans le fichier
        self.release_writer(RELEASE_FILE_PATH, species, release)

        # Configurer la barre de progression
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)

        # Lancer le téléchargement dans un thread séparé
        self.download_thread = DownloadThread(species, release)
        self.download_thread.progress.connect(self.progress_bar.setValue)
        self.download_thread.finished.connect(self.on_download_finished)
        self.download_thread.error.connect(self.on_download_error)
        self.download_thread.start()

    def on_download_finished(self):
        QMessageBox.information(self, "Download Complete", "Download and indexing complete.")
        self.progress_bar.setValue(100)
        self.accept()  # Ferme la boîte de dialogue

    def on_download_error(self, error):
        QMessageBox.critical(self, "Download Error", f"Failed to download Ensembl release: {error}")
        self.progress_bar.setVisible(False)

    def release_writer(self, file_path, species, release):
        with open(file_path, "w") as file:
            file.write(f"{species}\n{release}")

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
            ensembl_release.download()
            for i in range(1, 101, 10):  # Simulation de la progression
                self.msleep(200)  # Pause pour simuler un calcul
                self.progress.emit(i)
            ensembl_release.index()
            self.finished.emit()
        except Exception as e:
            self.error.emit(str(e))


def ask_ensembl_info():
    dialog = EnsemblDialog()
    if dialog.exec() == QDialog.DialogCode.Accepted:
        info = dialog.get_values()
        return info
    return None, None


if __name__ == "__main__":
    app = QApplication(sys.argv)
    species, release = ask_ensembl_info()
    if species and release:
        print(f"Espèce : {species}, Release : {release}")
    sys.exit(app.exec())
