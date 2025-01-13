WINDOW_HEIGHT = 600
WINDOW_WIDTH = 800

REF_COUPLE = ("start_genomic", "end_genomic")
A5SS_COL = ["shortSplice", "longSplice", "shareSplice"]
A3SS_COL = ["shortSplice", "longSplice", "shareSplice"]
RI_COL = ["RiStart", "RiEnd"]
SE_COL = ["upstreamEnd", "DownstreamStart"]
MXE_COL = ["1stExonStart", "1stExonEnd", "2ndExonStart", "2ndExonEnd", "upstreamEE", "downstreamES"]


import os 
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ICON_PATH = os.path.join(BASE_DIR, "Ressources", "Icones", "fugue-icons-3.5.6", "icons/")
QSS_PATH = os.path.join(BASE_DIR, "QSS", "styles.css")
RELEASE_FILE_PATH = os.path.join(BASE_DIR, "Ressources", "release", "release.txt")