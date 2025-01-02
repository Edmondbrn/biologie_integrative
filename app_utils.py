import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QWidget, QPushButton,
    QLabel, QToolBar, QStatusBar, QDialog, QFileDialog
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt



class FileDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Manual Calculation")
        self.setGeometry(100, 100, 400, 200)

        layout = QVBoxLayout()

        self.label = QLabel("Sélectionnez un fichier pour le calcul manuel des distances")
        layout.addWidget(self.label)

        self.button = QPushButton("Ouvrir l'explorateur de fichiers")
        self.button.clicked.connect(self.open_file_dialog)
        layout.addWidget(self.button)

        self.setLayout(layout)

   
    def open_file_dialog(self):
        file_dialog = QFileDialog()
        file_dialog.setNameFilter("Tableau CSV (*.csv);;Tableau TSV (*.tsv);;Tableau Excel (*.xlsx)")
        file_path, _ = file_dialog.getOpenFileName(self, "Open File", "", "Tableau CSV (*.csv);;Tableau TSV (*.tsv);;Tableau Excel (*.xlsx)")
        if file_path:
            self.label.setText(f"Fichier sélectionné : {file_path}")