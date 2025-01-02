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
        """
        Method that opens a file dialog to select a file.
        There are filters for CSV, TSV and Excel files.
        """
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "Open File", "", "CSV Table (*.csv);;TSV table (*.tsv);;Excel Table (*.xlsx)")
        if file_path:
            self.label.setText(f"Fichier sélectionné : {file_path}")