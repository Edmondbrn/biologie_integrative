from PyQt6.QtWidgets import QApplication
from src.Scripts.GUI.app import MainWindow
import sys
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())