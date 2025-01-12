from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel
from PyQt6.QtCore import Qt

class SimpleWindow(QWidget):
    def __init__(self, message):
        super().__init__()

        # Configurer la fenêtre
        self.setWindowTitle("Erreur de chargement")
        self.setGeometry(574, 400, 300, 100)

        # Créer un layout vertical
        layout = QVBoxLayout()

        # Ajouter un label avec le message
        label = QLabel(message)
        label.setStyleSheet("font-size: 16px; color: white;")  # Optionnel : style du texte
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)  # Centrer le texte

        # Ajouter le label au layout
        layout.addWidget(label)

        # Appliquer le layout à la fenêtre
        self.setLayout(layout)