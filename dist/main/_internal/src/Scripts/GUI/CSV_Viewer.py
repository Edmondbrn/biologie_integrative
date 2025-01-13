from PyQt6.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem
from PyQt6.QtGui import QIcon

from PyQt6.QtWidgets import (
    QMainWindow, QVBoxLayout, QWidget, QPushButton,
    QToolBar, QStatusBar, QMenu, QSpacerItem, QSizePolicy,
    QFileDialog, QToolButton, QTableWidgetItem, QTableWidget, 
)
from PyQt6.QtGui import QAction, QIcon
from PyQt6.QtCore import Qt, QSize, QPoint

from .manual_distances_window import ManualDistancesWindow
from .splicing_distances_window import SplicingDistancesWindow
from .all_splicing_distances_window import AllSplicingDistancesWindow
from .app_utils import show_alert, load_stylesheet

from ..Back.DrawGene import GeneImage
from ..GLOBAL import *

import pandas as pd
import pyensembl as pb

class CSVViewer(QWidget):
    def __init__(self, file_object, file_name):
        super().__init__()
        self.setWindowIcon(QIcon(f"{ICON_PATH}BI_logo.png"))
        self.setWindowTitle(file_name)
        self.setGeometry(100, 100, 800, 600)
        self.species, self.release = self.release_reader(RELEASE_FILE_PATH)
        self.release = int(self.release)

        layout = QVBoxLayout(self)
        self.tableWidget = QTableWidget()
        layout.addWidget(self.tableWidget)

        data = file_object.values.tolist()
        headers = file_object.columns.tolist()

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
        
        specy = self.species
        release = self.release
        try:
            bdd = pb.EnsemblRelease(species = specy, release = release)
            transcript : pb.Transcript = bdd.transcript_by_id(transcript_id)
            exon_pos = transcript.exon_intervals
            marker_pos = [min(texts[0], texts[1]), max(texts[0], texts[1])]
            gene = GeneImage(exon_pos, 
                            marker_pos, 
                            show_labels=show_distance_labels,
                            bar_xmax=marker_pos[1], 
                            bar_xmin=marker_pos[0])
            gene.show()
        except Exception as e:
            show_alert("Error", f"An error occured: {e}")
            return
        
    def release_reader(self, file_path):
        lines = []
        with open(file_path, "r") as file:
            for line in file:
                lines.append(line.strip())
            return lines