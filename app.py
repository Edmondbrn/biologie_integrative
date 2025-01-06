import sys
from PyQt6.QtWidgets import (
    QMainWindow, QApplication, QVBoxLayout, QWidget, QPushButton,
    QToolBar, QStatusBar, QMenu, QSpacerItem, QSizePolicy,
    QFileDialog, QToolButton, QTableWidgetItem, QTableWidget, 
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt, QSize, QPoint
from manual_distances_window import ManualDistancesWindow
from splicing_distances_window import SplicingDistancesWindow
from all_splicing_distances_window import AllSplicingDistancesWindow
from GLOBAL import *

from distances import *
from app_utils import load_stylesheet
from DrawGene import GeneImage
from app_utils import show_alert

ICON_PATH = "Ressources/Icones/fugue-icons-3.5.6/icons/"



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

        self.setStyleSheet(load_stylesheet("styles.css"))

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
        action_all_splicing.triggered.connect(lambda: self.onCalculateAllDistances())
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
        dialog = ManualDistancesWindow()
        dialog.exec()

    def onCalculateDistances(self, splice_type):
        dialog = SplicingDistancesWindow(splice_type)
        dialog.exec()

    def onCalculateAllDistances(self):
        dialog = AllSplicingDistancesWindow("all_splicing")
        dialog.exec()



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
            if verification_column in file_object.columns:
                self.prot_file = file_object
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
            print(self.genomic_file.head())
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

class CSVViewer(QWidget):
    def __init__(self, file_object):
        super().__init__()
        self.setWindowIcon(QIcon(f"{ICON_PATH}BI_logo.png"))
        self.setWindowTitle("CSV Viewer")
        self.setGeometry(100, 100, 800, 600)

        layout = QVBoxLayout(self)
        self.tableWidget = QTableWidget()
        layout.addWidget(self.tableWidget)

        data = file_object.values.tolist()
        headers = file_object.columns.tolist()
        if headers[0] == "Unnamed: 0":
            headers = headers[1:]

        # Remplir la table avec les données du CSV
        self.tableWidget.setRowCount(len(data))
        self.tableWidget.setColumnCount(len(data[0]))
        self.tableWidget.setHorizontalHeaderLabels(headers)

        for rowIdx, row in enumerate(data):
            for colIdx, cell in enumerate(row):
                self.tableWidget.setItem(rowIdx, colIdx, QTableWidgetItem(str(cell))) 
        self.tableWidget.resizeColumnsToContents()  # Adapter la taille des colonnes au contenu
        # -- AJOUT : autoriser le menu contextuel personnalisé --
        self.tableWidget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tableWidget.customContextMenuRequested.connect(self.onCustomContextMenu)

    def onCustomContextMenu(self, pos: QPoint):
        """
        This method is called when the user right-clicks on the QTableWidget.
        We retrieve the selected items and display a context menu.
        """
        selected_items = self.tableWidget.selectedItems()
        # definition du bouton pour appler la fonction
        menu = QMenu(self)
        action_graph_true = QAction("Display distances (labels ON)", self)
        action_graph_false = QAction("Display distances (lables OFF)", self)

        # Check if 3 cells are selected
        if len(selected_items) == 3:
            action_graph_true.setEnabled(True)
            action_graph_true.triggered.connect(lambda: self.show_graph(selected_items, show_distance_labels=True))
            action_graph_false.setEnabled(True)
            action_graph_false.triggered.connect(lambda: self.show_graph(selected_items, show_distance_labels=False))
        else:
            action_graph_true.setEnabled(False)
            action_graph_false.setEnabled(False)

        menu.addAction(action_graph_false)
        menu.addAction(action_graph_true)
        menu.exec(self.tableWidget.mapToGlobal(pos))

    def show_graph(self, items, transcript_id : str = None, show_distance_labels : bool = True):
        """
        Méthode appelée lorsque 3 cellules sont sélectionnées 
        et que l'utilisateur clique sur l'action dans le menu contextuel.
        """
        texts = [item.text() for item in items] # get the cell values
        for data in texts: # searching for the transcript ID in the selected cells
            if "ENS"  in data:
                transcript_id = data
                break
        if transcript_id is None:
            show_alert("Error", "Please select a transcript ID")
            return
        try:
            texts.remove(transcript_id) # remove the transcript ID from the list before the int casting
            texts = [int(text) for text in texts]
        except:
            show_alert("Error", "Please select cells with integer values except the transcript ID")
            return
        
        # TODO à changer pour les versions des génoms et les espèces
        specy = "mus_musculus"
        release = "102"
        bdd = pb.EnsemblRelease(species = specy, release = release)
        transcript : pb.Transcript = bdd.transcript_by_id(transcript_id)
        exon_pos = transcript.exon_intervals
        marker_pos = [min(texts[0], texts[1]), max(texts[0], texts[1])]
        gene = GeneImage(exon_pos, 
                         marker_pos, 
                         show_labels=show_distance_labels,
                         marker_colors=['red', 'blue'], 
                         bar_xmax=marker_pos[1], 
                         bar_xmin=marker_pos[0])
        gene.show()
        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())