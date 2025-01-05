
from manual_distances_window import ManualDistancesWindow
from PyQt6.QtWidgets import QLabel, QPlainTextEdit, QPushButton, QHBoxLayout
from app_utils import show_alert
from GLOBAL import *
class SplicingDistancesWindow(ManualDistancesWindow):
    
    def __init__(self, splice_type : str):
        super().__init__()
        self.splice = splice_type
        self.couple = []

    def addColumnsSelection(self):
        """
        Method to add the column selection widgets for comparing columns from the two dataframes.
        """
        self.column_selection_label = QLabel("Select columns to compare:")
        self.group_layout.addWidget(self.column_selection_label)

        # Zone de texte pour afficher les paires
        self.comparison_label = QLabel("Comparison pairs:")
        self.comparison_text = QPlainTextEdit()
        self.comparison_text.setReadOnly(True)
        self.group_layout.addWidget(self.comparison_text)

        self.button_compare_box = QHBoxLayout()
        # Bouton "Compare"
        self.compare_button = QPushButton("Compare")
        self.compare_button.clicked.connect(self.compare_columns)
        self.button_compare_box.addWidget(self.compare_button)

        # Fixer les couples de colonnes à comparer
        for couple in self.GenerateCouple():
            pair = f"{couple[0]} - {couple[1]}"
            self.comparison_text.appendPlainText(pair)
            self.compare_pairs.append(pair) # liste pour la classe mère pour qu'elle le passe au worker
        self.group_layout.addLayout(self.button_compare_box)

    def GenerateCouple(self):
        """
        Method to generate the couple of columns to compare.
        """
        switcher = {
            "A5SS": self.handle_A5SS,
            "RI": self.handle_RI,
            "A3SS": self.handle_A3SS,
            "SE": self.handle_SE,
            "MXE": self.handle_MXE}
        func = switcher.get(self.splice, lambda : show_alert("Error", "Invalid splicing type."))
        return func()
    
    def handle_A5SS(self):
        """
        Method to handle the A5SS+ splicing type.
        """
        A5SS_plus_col = A5SS_COL
        for col_name_ref in REF_COUPLE:
            for col_name in A5SS_plus_col:
                self.couple.append((col_name_ref, col_name))
        return self.couple
    
    def handle_A3SS(self):
        """
        Method to handle the A3SS+ splicing type.
        """
        A3SS_plus_col = A3SS_COL
        for col_name_ref in REF_COUPLE:
            for col_name in A3SS_plus_col:
                self.couple.append((col_name_ref, col_name))
        return self.couple
    
    def handle_RI(self):
        """
        Method to handle the RI splicing type.
        """
        RI_col = RI_COL
        for col_name_ref in REF_COUPLE:
            for col_name in RI_col:
                self.couple.append((col_name_ref, col_name))
        return self.couple
    
    def handle_SE(self):
        """
        Method to handle the SE splicing type.
        """
        SE_col = SE_COL
        for col_name_ref in REF_COUPLE:
            for col_name in SE_col:
                self.couple.append((col_name_ref, col_name))
        return self.couple
    
    def handle_MXE(self):
        """
        Method to handle the MXE splicing type.
        """
        MXE_col = MXE_COL
        for col_name_ref in REF_COUPLE:
            for col_name in MXE_col:
                self.couple.append((col_name_ref, col_name))
        return self.couple