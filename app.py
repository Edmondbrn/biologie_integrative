import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QWidget, QPushButton,
    QLabel, QToolBar, QStatusBar, QCheckBox, QMenu
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

 
        # ============================ définition des sous menus ===============================
        calculate_distances_menu = QMenu("Calculate distances", self)
        calculate_distances_menu.setIcon(QIcon(f"{ICON_PATH}ruler.png"))
        # Ajoute les actions pour les différents types d'épissage alternatif
        for splice_type in [ "A5SS", "A3SS", "RI", "MXE", "SE"]:
            action = QAction(f"Calculate distances for {splice_type}", self)
            action.setStatusTip(f"Calculate distances for {splice_type} events")
            action.triggered.connect(lambda checked, st=splice_type: self.onCalculateDistances(st))
            calculate_distances_menu.addAction(action)

        # Ajoute les options all_splicing et manual
        action_all_splicing = QAction("Calculate distances for all splicing", self)
        action_all_splicing.triggered.connect(lambda: self.onCalculateDistances("all_splicing"))
        action_all_splicing.setStatusTip("Calculate distances for all splicing events")
        calculate_distances_menu.addAction(action_all_splicing)

        action_manual = QAction("Manual calculation", self)
        action_manual.triggered.connect(lambda: self.onCalculateDistances("manual"))
        action_manual.setStatusTip("Calculate distances manually")
        calculate_distances_menu.addAction(action_manual)


        self.setStatusBar(QStatusBar(self))

        menu = self.menuBar()
        # menu File de la barre de tâches
        file_menu = menu.addMenu("File")
        file_menu.addAction(button_load_data)
        file_menu.addAction(buttton_save_data)
        file_menu.addAction(button_quit)
        # Menu Action de la barre de tâches
        action_menu = menu.addMenu("Actions")
        action_menu.addAction(button_convert)
        action_menu.addMenu(calculate_distances_menu)





    def onMyToolBarButtonClick(self, s):
        print("Button clicked", s)
    def onCalculateDistances(self, splice_type):
        print(f"Calculating distances for {splice_type}")


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()