from PyQt6.QtWidgets import QMessageBox

def load_stylesheet(file_path):
    with open(file_path, "r") as file:
        return file.read()
    

def show_alert(title: str, message: str):
    """
    Method to display an alert dialog with a title and a specific message.
    """
    dict_logo = {"Info": QMessageBox.Icon.Information, "Error": QMessageBox.Icon.Critical, "Warning": QMessageBox.Icon.Warning}
    alert = QMessageBox()
    alert.setWindowTitle(title)
    alert.setText(message)
    alert.setIcon(dict_logo[title])
    alert.setStandardButtons(QMessageBox.StandardButton.Ok)
    alert.exec()