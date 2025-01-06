from GLOBAL import ICON_PATH

from PyQt6.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem
from PyQt6.QtGui import QIcon

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
        self.tableWidget.setColumnCount(len(data[0])-1)
        self.tableWidget.setHorizontalHeaderLabels(headers)

        for rowIdx, row in enumerate(data):
            for colIdx, cell in enumerate(row[1:]):
                self.tableWidget.setItem(rowIdx, colIdx, QTableWidgetItem(str(cell))) 
