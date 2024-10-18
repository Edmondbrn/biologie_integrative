import pandas as pd
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
print(os.getcwd())

df = pd.read_csv("rmats_post/A3SS.MATS.JCEC.txt", sep="\t")