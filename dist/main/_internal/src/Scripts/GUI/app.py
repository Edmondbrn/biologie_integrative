from PyQt6.QtWidgets import (
    QMainWindow, QVBoxLayout, QWidget, QPushButton,
    QToolBar, QStatusBar, QMenu, QSpacerItem, QSizePolicy,
    QFileDialog, QToolButton, QTableWidgetItem, QTableWidget, 
    QApplication, QLabel, QWidgetAction
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt, QSize, QPoint

from .manual_distances_window import ManualDistancesWindow
from .splicing_distances_window import SplicingDistancesWindow
from .all_splicing_distances_window import AllSplicingDistancesWindow
from .app_utils import show_alert, load_stylesheet

from ..Back.DrawGene import GeneImage
from ..Back.Id_convertor import *
from .. import GLOBAL

import sys
import pandas as pd
import pyensembl as pb
import csv
import os

from .ID_convertor import IDConversionDialog
from .CSV_Viewer import CSVViewer

from .RNAtoDNA import RNAtoDNAWindow
from .parsingRmats import ParsingRmats
from ..GLOBAL import *
from .EnsemblDialog import EnsemblDialog

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()

        self.file_path = None
        self.prot_file = None
        self.genomic_file = None
        self.output_file = None 
        self.species, self.release = self.release_reader(GLOBAL.RELEASE_FILE_PATH)
        self.release = int(self.release)

        self.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)

        self.setWindowTitle("RepositionX")
        self.setWindowIcon(QIcon(f"{GLOBAL.ICON_PATH}BI_logo.png"))

        # Configure the window to start maximized
        # self.setWindowState(self.windowState() | Qt.WindowState.WindowMaximized)

        self.setStyleSheet(load_stylesheet(GLOBAL.QSS_PATH))

        self.__create_menu()
        # Central widget
        central_widget = QWidget()
        central_widget.setObjectName("central_widget")
        self.setCentralWidget(central_widget)

        # Set layout to central widget
        layout = QVBoxLayout(central_widget)
        central_widget.setLayout(layout)

        logo = QPushButton()
        logo.setIcon(QIcon(f"{GLOBAL.ICON_PATH}BI_logo.png"))
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
        button_load_data = QAction(QIcon(f"{GLOBAL.ICON_PATH}document-excel.png"), "Load data", self)
        button_load_data.setStatusTip("Load data from a CSV table")
        button_load_data.triggered.connect(self.open_file_dialog)
        toolbar.addAction(button_load_data)

        # bouton pour charger un output
        button_load_output = QAction(QIcon(f"{GLOBAL.ICON_PATH}document-excel-table.png"), "Load output", self)
        button_load_output.setStatusTip("Load output from a CSV table")
        button_load_output.triggered.connect(self.open_output_dialog)
        toolbar.addAction(button_load_output)

        # bouton pour quitter l'application
        button_quit = QAction(QIcon(f"{GLOBAL.ICON_PATH}door-open-out.png"), "Quit", self)
        button_quit.setStatusTip("Close application")
        button_quit.triggered.connect(self.close)
        toolbar.addAction(button_quit)

        # Bouton pour convertir des coordonnées ARNm en coordonnées génomiques
        button_convert = QAction(QIcon(f"{GLOBAL.ICON_PATH}arrow-circle.png"), "Convert mRNA to DNA", self)
        button_convert.setStatusTip("Convert RNA coordinates to genomic coordinates")
        button_convert.triggered.connect(lambda: self.onRNAtoDNA())

        # Bouton pour traiter les fichiers rMATS
        button_convert2 = QAction(QIcon(f"{ICON_PATH}magnifier.png"), "Parsing rMATS", self)
        button_convert2.setStatusTip("Open the rMATS Parser tool")
        button_convert2.triggered.connect(self.openParsingRmatsWindow)



        #Bouton pour convertir les ID d'ensembl à NCBI
        button_convert_ID = QAction(QIcon(f"{GLOBAL.ICON_PATH}address-book-blue.png"), "Change NCBI IDs to Ensembl IDs", self)
        button_convert_ID.setStatusTip("Change NCBI IDs to Ensembl IDs")
        button_convert_ID.triggered.connect(lambda: self.ID_loader())

        button_change_release = QAction(QIcon(f"{GLOBAL.ICON_PATH}arrow-circle-double.png"), "Change release and species", self)
        button_change_release.setStatusTip("Change the release and species of the ensembl genome reference")
        button_change_release.triggered.connect(lambda: self.change_release())
        

        # ============================ définition des sous menus ===============================
        calculate_distances_menu = QMenu("Calculate distances", self)
        calculate_distances_menu.setIcon(QIcon(f"{GLOBAL.ICON_PATH}ruler.png"))
        # Ajoute les actions pour les différents types d'épissage alternatif
        for splice_type in [ "A5SS", "A3SS", "RI", "MXE", "SE"]:
            action = QAction(f"Calculate distances for {splice_type}", self)
            action.setStatusTip(f"Calculate distances for {splice_type} events")
            action.triggered.connect(lambda checked, st=splice_type: self.onCalculateDistances(st, self.prot_file, self.genomic_file))
            calculate_distances_menu.addAction(action)

        # Ajoute les options all_splicing et manual
        action_all_splicing = QAction("Calculate distances for all splicing", self)
        action_all_splicing.triggered.connect(lambda: self.onCalculateAllDistances("all_splicing", self.prot_file, self.genomic_file))
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
        file_menu.addAction(button_load_output)
        file_menu.addAction(button_quit)
        # Menu Action de la barre de tâches
        action_menu = menu.addMenu("Actions")
        action_menu.addAction(button_convert_ID)
        action_menu.addAction(button_convert)
        action_menu.addAction(button_convert2)
        action_menu.addMenu(calculate_distances_menu)

        #menu download
        download_menu = menu.addMenu("Release")
        download_menu.addAction(button_change_release)

        self.release_menu = QMenu(self)
        self.menuBar().addMenu(self.release_menu)

        # Initialiser le texte du menu
        self.update_release_menu()

    def update_release_menu(self):
        # Mettre à jour le titre du menu avec les valeurs actuelles de GLOBAL
        self.release_menu.setTitle(f"Release: {self.release}, Species: {self.species}")

    def change_release(self):
        dialog = EnsemblDialog()
        if dialog.exec():  
            self.species, self.release = self.release_reader(GLOBAL.RELEASE_FILE_PATH)
            self.update_release_menu()

    def onManualDistances(self, reference_file=None, genomic_file=None):
        dialog = ManualDistancesWindow(reference_file, genomic_file)
        dialog.data_signal.connect(self.update_files)  # Connecter le signal au slot
        dialog.name_signal.connect(self.update_file_names)  # Connecter le signal des noms au slot
        dialog.exec()

    def update_files(self, files):
        self.prot_file, self.genomic_file = files
        self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def update_file_names(self, names):
        self.prot_file_name, self.genomic_file_name = names
        self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def open_file_dialog(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "File Explorer", "", "All Files (*);;CSV Files (*.csv)")
        if file_path:
            self.file_path = file_path
            self.file_loader()

    def ID_loader(self):
        dialog = IDConversionDialog()
        dialog.exec()

    def open_output_dialog(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "File Explorer", "", "CSV Files (*.csv)")
        if file_path:
            self.file_path_output = file_path
            self.file_loader_output()
        

    def file_loader_output(self):
        if self.file_path_output:
            if self.output_file is not None:
                show_alert("Error", "A genomic file has already been loaded \n Please delete it before loading another file")
            else:
                if self.file_path_output.endswith(".csv"):
                    sep=self.detect_separator(self.file_path_output)
                    index_column=self.index_column_detector(self.file_path_output, sep)
                    self.output_file = pd.read_csv(self.file_path_output, sep=sep, index_col=index_column)
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))
                    self.output_name = os.path.basename(self.file_path_output)
                else:
                    show_alert("Error", "The file format is not supported")
                    return

    def index_column_detector(self, file_path, sep):
        file = pd.read_csv(file_path, sep=sep, nrows=10)
        if file.columns[0] == "" or 'Unnamed: 0' in file.columns[0]:
            return 0
        else:
            return None
        
    def index_column_detector_excel(self, file_path):
        file = pd.read_excel(file_path, nrows=10)
        if file.columns[0] == "" or 'Unnamed: 0' in file.columns[0]:
            return 0
        else:
            return None

    def file_loader(self):
        if self.file_path:
            if self.file_path.endswith(".csv"):
                sep=self.detect_separator(self.file_path)
                index_column=self.index_column_detector(self.file_path, sep)
                file_object = pd.read_csv(self.file_path, sep=sep, index_col=index_column)
            elif self.file_path.endswith(".tsv"):
                file_object = pd.read_csv(self.file_path, sep="\t")
            elif self.file_path.endswith(".xlsx"):
                index_column = self.index_column_detector_excel(self.file_path)
                file_object = pd.read_excel(self.file_path, index_col=index_column)
            else:
                show_alert("Error", "The file format is not supported")
                return
            verification_column = 'GeneID'
            if verification_column == file_object.columns[0]:
                if self.genomic_file is not None:
                    show_alert("Error", "A genomic file has already been loaded \n Please delete it before loading another file")
                else:
                    self.genomic_file = file_object
                    self.genomic_file_name = os.path.basename(self.file_path)
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))
            else:
                if self.prot_file is not None:
                    show_alert("Error", "A reference file has already been loaded \n Please delete it before loading another file")
                else:
                    self.prot_file = file_object
                    self.prot_file_name = os.path.basename(self.file_path)
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def close_file_custom(self, variable_name: str):
        setattr(self, variable_name, None)
        self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def release_reader(self, file_path):
        lines = []
        with open(file_path, "r") as file:
            for line in file:
                lines.append(line.strip())
            return lines

    def detect_separator(self, file_path):
        """
        Detect the separator used in a CSV file.
        """
        with open(file_path, 'r') as file:
        # Read the first line of the file
            sample = file.read(1024)
            # Try different separators
            sniffer = csv.Sniffer()
            try:
                sep = sniffer.sniff(sample).delimiter
            except csv.Error:
                # Default to tab if no separator is found
                sep = '\t'
            return sep
    
    def onCalculateDistances(self, splice_type, reference_file, genomic_file):
        dialog = SplicingDistancesWindow(splice_type, reference_file, genomic_file)
        dialog.data_signal.connect(self.update_files)  # Connecter le signal au slot
        dialog.name_signal.connect(self.update_file_names)  # Connecter le signal des noms au slot
        dialog.exec()

    def onCalculateAllDistances(self, splicing, reference_file, genomic_file):
        dialog = AllSplicingDistancesWindow(splicing, reference_file, genomic_file)
        dialog.data_signal.connect(self.update_files)  # Connecter le signal au slot
        dialog.name_signal.connect(self.update_file_names)  # Connecter le signal des noms au slot
        dialog.exec()

    def dynamic_menues(self, toolbar):
        if self.prot_file is not None and self.findChild(QMenu, "Fichier protéine") is None: 
            tool_button = QToolButton(self)
            tool_button.setText("reference file  ")
            tool_button.setObjectName("toolButtonProtein")
            tool_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)  # Le menu apparaît immédiatement

            # Créer un menu pour le QToolButton
            menu_protein = QMenu("Menu protéine", self)
            menu_protein.setObjectName("Fichier protéine")
            view_protein_button = QAction("View file", self)
            view_protein_button.setStatusTip("View file button")
            view_protein_button.triggered.connect(lambda: self.csv_viewer(self.prot_file, self.prot_file_name))
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
            view_genomic_button.triggered.connect(lambda: self.csv_viewer(self.genomic_file, self.genomic_file_name))
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

        if self.output_file is not None and self.findChild(QMenu, "output") is None: 
            tool_button_output = QToolButton(self)
            tool_button_output.setText("output file  ")
            tool_button_output.setObjectName("toolButtonOutput")
            tool_button_output.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)  # Le menu apparaît immédiatement

            # Créer un menu pour le QToolButton
            menu_output = QMenu("output", self)
            menu_output.setObjectName("output")
            view_output_button = QAction("View file", self)
            view_output_button.setStatusTip("View file button")
            view_output_button.triggered.connect(lambda: self.csv_viewer(self.output_file, self.output_name))
            menu_output.addAction(view_output_button)

            # Ajouter des actions au menu
            delete_output = QAction("delete", self)
            delete_output.setStatusTip("delete")
            delete_output.triggered.connect(lambda: self.close_file_custom("output_file"))
            menu_output.addAction(delete_output)

            # Associer le menu au QToolButton
            tool_button_output.setMenu(menu_output)

            # Ajouter le QToolButton à la barre d'outils
            toolbar.addWidget(tool_button_output)

        elif self.genomic_file is None and self.findChild(QMenu, "Fichier génomique") != None:
            self.findChild(QMenu, "Fichier génomique").deleteLater()
            self.findChild(QToolButton, "toolButtonGenomic").deleteLater()
        
        elif self.prot_file is None and self.findChild(QMenu, "Fichier protéine") != None:
            self.findChild(QMenu, "Fichier protéine").deleteLater()
            self.findChild(QToolButton, "toolButtonProtein").deleteLater()

        elif self.output_file is None and self.findChild(QMenu, "output") != None:
            self.findChild(QMenu, "output").deleteLater()
            self.findChild(QToolButton, "toolButtonOutput").deleteLater()

    def csv_viewer(self, file, file_name="CSV Viewer"):
        self.viewer = CSVViewer(file, file_name)
        self.viewer.show()

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
