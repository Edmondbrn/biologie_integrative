import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QPushButton, QLabel, QVBoxLayout
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap

def load_stylesheet(file_path):
    with open(file_path, "r") as file:
        return file.read()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # Configure the window to start maximized
        self.setWindowState(self.windowState() | Qt.WindowMaximized)

        # Central widget
        central_widget = QWidget()
        central_widget.setObjectName("central_widget")
        self.setCentralWidget(central_widget)

        self.setStyleSheet(load_stylesheet("styles.qss"))

    def add_non_central_widget(self):
        # Create a widget that will not be central
        non_central_widget = QWidget(self)

        # Layout for the non-central widget
        non_central_layout = QVBoxLayout(non_central_widget)

        # Create a button for demonstration
        button = QPushButton("Click me", non_central_widget)
        non_central_layout.addWidget(button)

        # Set the non-central widget's layout
        non_central_widget.setLayout(non_central_layout)

        # Position the widget in the main window
        non_central_widget.move(100, 100)  # You can adjust these coordinates for responsiveness

        # Resize the widget to take 20% of the window's width and height
        non_central_widget.resize(self.width() * 0.2, self.height() * 0.2)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
