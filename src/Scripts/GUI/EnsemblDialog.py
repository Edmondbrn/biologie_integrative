from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QApplication

class EnsemblDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Sélection de l'espèce et release Ensembl")

        layout = QVBoxLayout(self)

        # Champ pour l'espèce
        species_layout = QHBoxLayout()
        species_label = QLabel("Specy :")
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
        btn_ok.clicked.connect(self.accept)
        layout.addWidget(btn_ok)

    def get_values(self):
        return self.species_input.text(), self.release_input.text()
    
    
def ask_ensembl_info():
    dialog = EnsemblDialog()
    if dialog.exec() == QDialog.DialogCode.Accepted:
        return dialog.get_values()
    return None, None

if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    species, release = ask_ensembl_info()
    if species and release:
        print(f"Espèce : {species}, Release : {release}")
    sys.exit(app.exec())