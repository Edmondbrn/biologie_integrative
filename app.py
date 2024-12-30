import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QWidget, QPushButton,
    QLabel, QToolBar, QStatusBar, QCheckBox
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt

ICON_PATH = "Ressources/Icones/fugue-icons-3.5.6/icons/"

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Biologie intégrative")

        layout = QVBoxLayout()
        button = QPushButton("Cliquez-moi")
        layout.addWidget(button)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)
        self.__create_menu()

        self.resize(800, 600)


    def __create_menu(self):
        # On créé une barre d'outils
        toolbar = QToolBar("My main toolbar")
        self.addToolBar(toolbar)

        # ============================ définition des boutons ===============================
        # bouton pour charger un fichier
        button_load_data = QAction(QIcon(f"{ICON_PATH}document-excel.png"), "Load data", self)
        button_load_data.setStatusTip("Load data from a CSV table")
        button_load_data.triggered.connect(lambda: self.onMyToolBarButtonClick("test pression bouton"))
        toolbar.addAction(button_load_data)

        # Bouton pour sauvegarder les données
        buttton_save_data = QAction(QIcon(f"{ICON_PATH}disk.png"), "Save data", self)
        buttton_save_data.setStatusTip("Save data to a CSV table")
        buttton_save_data.triggered.connect(lambda: self.onMyToolBarButtonClick("test pression bouton"))
        toolbar.addAction(buttton_save_data)

        # bouton pour quitter l'application
        button_quit = QAction(QIcon(f"{ICON_PATH}door-open-out.png"), "Quit", self)
        button_quit.setStatusTip("Close application")
        button_quit.triggered.connect(self.close)
        toolbar.addAction(button_quit)

        # Bouton pour convertir des coordonnées ARNm en coordonnées génomiques
        button_convert = QAction(QIcon(f"{ICON_PATH}arrow-circle.png"), "Convert mRNA to DNA", self)
        button_convert.setStatusTip("Convert RNA coordinates to genomic coordinates")
        button_convert.triggered.connect(lambda: self.onMyToolBarButtonClick("test pression bouton"))

        # Bouton pour calculer les distances entre les coordonnées
        button_distances = QAction(QIcon(f"{ICON_PATH}ruler.png"), "Calculate distances", self)
        button_distances.setStatusTip("Calculate distances between coordinates")
        button_distances.triggered.connect(lambda: self.onMyToolBarButtonClick("test pression bouton"))



        self.setStatusBar(QStatusBar(self))

        menu = self.menuBar()
        file_menu = menu.addMenu("File")
        file_menu.addAction(button_load_data)
        file_menu.addAction(buttton_save_data)
        file_menu.addAction(button_quit)

        action_menu = menu.addMenu("Actions")
        action_menu.addAction(button_convert)
        action_menu.addAction(button_distances)




    def onMyToolBarButtonClick(self, s):
        print("Button clicked", s)

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()