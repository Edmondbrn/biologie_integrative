from PyQt6.QtWidgets import (
    QMainWindow, QVBoxLayout, QWidget, QPushButton,
    QToolBar, QStatusBar, QMenu, QSpacerItem, QSizePolicy,
    QFileDialog, QToolButton, QTableWidgetItem, QTableWidget, QApplication
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt, QSize, QPoint

from .manual_distances_window import ManualDistancesWindow
from .splicing_distances_window import SplicingDistancesWindow
from .all_splicing_distances_window import AllSplicingDistancesWindow
from .app_utils import show_alert, load_stylesheet

from ..Back.DrawGene import GeneImage
from ..GLOBAL import *

import sys
import pandas as pd
import pyensembl as pb

from ..Back.distances import *
from .CSV_Viewer import CSVViewer

from .RNAtoDNA import RNAtoDNAWindow
from .parsingRmats import ParsingRmats


def load_stylesheet(file_path):
    with open(file_path, "r") as file:
        return file.read()

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()

        self.file_path = None
        self.prot_file = None
        self.genomic_file = None

        self.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)

        self.setWindowTitle("BI Project")
        self.setWindowIcon(QIcon(f"{ICON_PATH}BI_logo.png"))

        # Configure the window to start maximized
        self.setWindowState(self.windowState() | Qt.WindowState.WindowMaximized)

        self.setStyleSheet(load_stylesheet(QSS_PATH))

        self.__create_menu()

        # Central widget
        central_widget = QWidget()
        central_widget.setObjectName("central_widget")
        self.setCentralWidget(central_widget)

        # Set layout to central widget
        layout = QVBoxLayout(central_widget)
        central_widget.setLayout(layout)

        logo = QPushButton()
        logo.setIcon(QIcon(f"{ICON_PATH}BI_logo.png"))
        logo.setObjectName('logo')
        logo.setIconSize(QSize(500, 500))  # Set the icon size
        logo.setFixedSize(500, 500)
        layout.addWidget(logo)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)

        spacer = QSpacerItem(20, 100, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
        layout.addItem(spacer)

    def __create_menu(self):
        # On créé une barre d'outils
        toolbar = QToolBar("My main toolbar")
        toolbar.setObjectName("My main toolbar")
        self.addToolBar(toolbar)

        # ============================ définition des boutons ===============================
        # bouton pour charger un fichier
        button_load_data = QAction(QIcon(f"{ICON_PATH}document-excel.png"), "Load data", self)
        button_load_data.setStatusTip("Load data from a CSV table")
        button_load_data.triggered.connect(self.open_file_dialog)
        toolbar.addAction(button_load_data)

        # bouton pour quitter l'application
        button_quit = QAction(QIcon(f"{ICON_PATH}door-open-out.png"), "Quit", self)
        button_quit.setStatusTip("Close application")
        button_quit.triggered.connect(self.close)
        toolbar.addAction(button_quit)

        # Bouton pour convertir des coordonnées ARNm en coordonnées génomiques
        button_convert = QAction(QIcon(f"{ICON_PATH}arrow-circle.png"), "Convert mRNA to DNA", self)
        button_convert.setStatusTip("Convert RNA coordinates to genomic coordinates")
        button_convert.triggered.connect(lambda: self.onRNAtoDNA())

        # Bouton pour traiter les fichiers rMATS
        button_convert2 = QAction(QIcon(f"{ICON_PATH}magnifier.png"), "Parsing rMATS", self)
        button_convert2.setStatusTip("Open the rMATS Parser tool")
        button_convert2.triggered.connect(self.openParsingRmatsWindow)



 
        # ============================ définition des sous menus ===============================
        calculate_distances_menu = QMenu("Calculate distances", self)
        calculate_distances_menu.setIcon(QIcon(f"{ICON_PATH}ruler.png"))
        # Ajoute les actions pour les différents types d'épissage alternatif
        for splice_type in [ "A5SS", "A3SS", "RI", "MXE", "SE"]:
            action = QAction(f"Calculate distances for {splice_type}", self)
            action.setStatusTip(f"Calculate distances for {splice_type} events")
            action.triggered.connect(lambda checked, st=splice_type: self.onCalculateDistances(st, self.prot_file, self.genomic_file))
            calculate_distances_menu.addAction(action)

        # Ajoute les options all_splicing et manual
        action_all_splicing = QAction("Calculate distances for all splicing", self)
        action_all_splicing.triggered.connect(lambda: self.onCalculateDistances("all_splicing", self.prot_file, self.genomic_file))
        action_all_splicing.setStatusTip("Calculate distances for all splicing events")
        calculate_distances_menu.addAction(action_all_splicing)

        action_manual = QAction("Manual calculation", self)
        action_manual.triggered.connect(lambda: self.onManualDistances(self.prot_file, self.genomic_file))
        action_manual.setStatusTip("Calculate distances manually")
        calculate_distances_menu.addAction(action_manual)


        self.setStatusBar(QStatusBar(self))

        menu = self.menuBar()
        # menu File de la barre de tâches
        file_menu = menu.addMenu("File")
        file_menu.addAction(button_load_data)
        file_menu.addAction(button_quit)
        # Menu Action de la barre de tâches
        action_menu = menu.addMenu("Actions")
        action_menu.addAction(button_convert)
        action_menu.addAction(button_convert2)
        action_menu.addMenu(calculate_distances_menu)

    def onManualDistances(self, reference_file=None, genomic_file=None):
        dialog = ManualDistancesWindow(reference_file, genomic_file)
        dialog.data_signal.connect(self.update_files)  # Connecter le signal au slot
        dialog.exec()

    def update_files(self, files):
        self.prot_file, self.genomic_file = files
        self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def open_file_dialog(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "File Explorer", "", "All Files (*);;CSV Files (*.csv)")
        if file_path:
            self.file_path = file_path
            self.file_loader()

    def save_file_dialog(self):
        file_name, _ = QFileDialog.getSaveFileName(self, "File Explorer", "", "All Files (*);;CSV Files (*.csv)")
        if file_name:
            self.file_name = file_name

    def file_loader(self):
        if self.file_path:
            if self.file_path.endswith(".csv"):
                file_object = pd.read_csv(self.file_path, sep=self.detect_separator(self.file_path), index_col=0)
            elif self.file_path.endswith(".tsv"):
                file_object = pd.read_csv(self.file_path, sep="\t")
            elif self.file_path.endswith(".xlsx"):
                file_object = pd.read_excel(self.file_path, index_col=0)
            verification_column = 'GeneID'
            if verification_column == file_object.columns[0]: #TODO vérification sur le type de fichier pour pas ouvrir n'importe quoi, et faire un tri également sur les colonnes du fichier génomique
                if self.genomic_file is not None:
                    self.error_window("A genomic file has already been loaded\n Please delete it before loading another file") 
                else:
                    self.genomic_file = file_object
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))
            else:
                if self.prot_file is not None:
                    self.error_window("A reference file has already been loaded \n Please delete it before loading another file")
                else:
                    self.prot_file = file_object
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def close_file_custom(self, variable_name: str):
        setattr(self, variable_name, None)
        self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))


    def detect_separator(self, file_path):
        """
        Detect the separator used in a CSV file.
        """
        with open(file_path, 'r') as file:
            # Read the first line of the file
            first_line = file.readline()
            # Try different separators
            separators = [',', ';', '\t', ' ']
            for sep in separators:
                if sep in first_line:
                    return sep
        return ','  # Default to comma if no separator is found
    
    def onCalculateDistances(self, splice_type, reference_file, genomic_file):
        dialog = SplicingDistancesWindow(splice_type, reference_file, genomic_file)
        dialog.exec()

    def onCalculateAllDistances(self, reference_file, genomic_file):
        dialog = AllSplicingDistancesWindow("all_splicing", reference_file, genomic_file)
        dialog.exec()

    def dynamic_menues(self, toolbar):
        if self.prot_file is not None and self.findChild(QMenu, "Fichier protéine") is None: #TODO n'arrive pas à lire les fichiers pour FMRP
            tool_button = QToolButton(self)
            tool_button.setText("reference file  ")
            tool_button.setObjectName("toolButtonProtein")
            tool_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)  # Le menu apparaît immédiatement

            # Créer un menu pour le QToolButton
            menu_protein = QMenu("Menu protéine", self)
            menu_protein.setObjectName("Fichier protéine")
            view_protein_button = QAction("View file", self)
            view_protein_button.setStatusTip("View file button")
            view_protein_button.triggered.connect(lambda: self.csv_viewer(self.prot_file))
            menu_protein.addAction(view_protein_button)

            # Ajouter des actions au menu
            delete_protein = QAction("delete", self)
            delete_protein.setStatusTip("delete")
            delete_protein.triggered.connect(lambda: self.close_file_custom("prot_file"))
            menu_protein.addAction(delete_protein)

            # Associer le menu au QToolButton
            tool_button.setMenu(menu_protein)

            # Ajouter le QToolButton à la barre d'outils
            toolbar.addWidget(tool_button)

        if self.genomic_file is not None and self.findChild(QMenu, "Fichier génomique") is None:
            tool_button_genomic = QToolButton(self)
            tool_button_genomic.setText("genomic file  ")
            tool_button_genomic.setObjectName("toolButtonGenomic")
            tool_button_genomic.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)  # Le menu apparaît immédiatement

            # Créer un menu pour le QToolButton
            menu_genomic = QMenu("Fichier génomique", self)
            menu_genomic.setObjectName("Fichier génomique")
            view_genomic_button = QAction("View file", self)
            view_genomic_button.setStatusTip("View file button")
            view_genomic_button.triggered.connect(lambda: self.csv_viewer(self.genomic_file))
            menu_genomic.addAction(view_genomic_button)

            # Ajouter des actions au menu
            delete_genomic = QAction("delete", self)
            delete_genomic.setStatusTip("delete")
            delete_genomic.triggered.connect(lambda: self.close_file_custom("genomic_file"))
            menu_genomic.addAction(delete_genomic)

            # Associer le menu au QToolButton
            tool_button_genomic.setMenu(menu_genomic)

            # Ajouter le QToolButton à la barre d'outils
            toolbar.addWidget(tool_button_genomic)

        elif self.genomic_file is None and self.findChild(QMenu, "Fichier génomique") != None:
            self.findChild(QToolButton, "toolButtonGenomic").deleteLater()
        
        elif self.prot_file is None and self.findChild(QMenu, "Fichier protéine") != None:
            self.findChild(QToolButton, "toolButtonProtein").deleteLater()

    def csv_viewer(self, file_name):
        self.viewer = CSVViewer(file_name)
        self.viewer.show()
    
    def error_window(self, message):
        self.error = SimpleWindow(message)
        self.error.show()
    
    def onRNAtoDNA(self):
        """
        Cette méthode est appelée lorsque l'on clique sur le bouton dans la barre d'outils.
        Elle ouvrira la fenêtre RNAtoDNAWindow pour convertir les coordonnées mRNA en coordonnées génomiques.
        """
        # Créez l'instance de la fenêtre RNAtoDNAWindow
        window = RNAtoDNAWindow()
        window.exec()  # Affichez la fenêtre

    def openParsingRmatsWindow(self):
    # Création et affichage de la fenêtre ParsingRmats
        self.parsing_window = ParsingRmats()  # Instancie la fenêtre ParsingRmats
        self.parsing_window.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
