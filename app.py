import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QWidget, QPushButton,
    QLabel, QToolBar, QStatusBar, QCheckBox, QMenu, QSpacerItem, QSizePolicy,
    QFileDialog, QToolButton, QTableWidgetItem, QTableWidget, 
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt, QSize
from app_utils import FileDialogManual
from GLOBAL import *

from distances import *
from CSV_Viewer import CSVViewer
from Error_window import SimpleWindow

def load_stylesheet(file_path):
    with open(file_path, "r") as file:
        return file.read()

class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()

        self.file_path = None
        self.prot_file = None
        self.genomic_file = None

        self.setWindowTitle("BI Project")
        self.setWindowIcon(QIcon(f"{ICON_PATH}BI_logo.png"))

        # Configure the window to start maximized
        self.setWindowState(self.windowState() | Qt.WindowState.WindowMaximized)

        self.setStyleSheet(load_stylesheet("styles.qss"))

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
        action_manual.triggered.connect(lambda: self.onManualDistances())
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

    def onManualDistances(self):
        dialog = FileDialogManual()
        dialog.exec()


    def onMyToolBarButtonClick(self, s):
        print("Button clicked", s)
    def onCalculateDistances(self, splice_type):
        print(f"Calculating distances for {splice_type}")

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
            file_object = pd.read_csv(self.file_path, sep="\t")
            verification_column = 'gene_name'
            if verification_column in file_object.columns: #TODO vérification sur le type de fichier pour pas ouvrir n'importe quoi, et faire un tri également sur les colonnes du fichier génomique
                if self.prot_file is not None:
                    self.error_window("Un fichier protéine a déjà été chargé. \n Veuillez le supprimer avant de charger un nouveau fichier")
                else:
                    self.prot_file = file_object
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))
            else:
                if self.genomic_file is not None:
                    self.error_window("Un fichier génomique a déjà été chargé. \n Veuillez le supprimer avant de charger un nouveau fichier") 
                else:
                    self.genomic_file = file_object
                    self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))

    def close_file_custom(self, variable_name: str):
        setattr(self, variable_name, None)
        self.dynamic_menues(self.findChild(QToolBar, "My main toolbar"))


    def dynamic_menues(self, toolbar):
        if self.prot_file is not None: #TODO n'arrive pas à lire les fichiers pour FMRP
            tool_button = QToolButton(self)
            tool_button.setText("Fichier protéine  ")
            tool_button.setObjectName("toolButtonProtein")
            tool_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)  # Le menu apparaît immédiatement

            # Créer un menu pour le QToolButton
            menu_protein = QMenu("Menu protéine", self)

            # Ajouter des actions au menu
            delete_prot = QAction("Supprimer fichier protéine", self)
            delete_prot.setStatusTip("Supprimer le fichier protéine")
            delete_prot.triggered.connect(lambda: self.close_file_custom("prot_file"))
            menu_protein.addAction(delete_prot)

            # Associer le menu au QToolButton
            tool_button.setMenu(menu_protein)

            # Ajouter le QToolButton à la barre d'outils
            toolbar.addWidget(tool_button)

        if self.genomic_file is not None:
            tool_button_genomic = QToolButton(self)
            tool_button_genomic.setText("Fichier génomique  ")
            tool_button_genomic.setObjectName("toolButtonGenomic")
            tool_button_genomic.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)  # Le menu apparaît immédiatement

            # Créer un menu pour le QToolButton
            menu_genomic = QMenu("Fichier génomique", self)
            menu_genomic.setObjectName("Fichier génomique")
            view_genomic_button = QAction("Visualiser le fichier", self)
            view_genomic_button.setStatusTip("Visualiser le fichier génomique")
            view_genomic_button.triggered.connect(lambda: self.csv_viewer(self.genomic_file))
            menu_genomic.addAction(view_genomic_button)

            # Ajouter des actions au menu
            delete_genomic = QAction("Supprimer fichier génomique", self)
            delete_genomic.setStatusTip("Supprimer le fichier génomique")
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


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())