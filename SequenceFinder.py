import pyensembl as pb
from pandas import read_csv
import os
import re
from multiprocessing import Pool
import time


os.chdir(os.path.dirname(os.path.abspath(__file__)))

class SequenceFinder():
    """
    Cette classe va permettre de récupérer
    """
