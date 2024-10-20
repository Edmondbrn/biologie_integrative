import pandas as pd
from GeneFinder import GeneFinder
import numpy as np
import os

os.chdir(os.path.dirname(__file__))


fh = open("FMRP_Binding_sites_mouse_Maurin_NAR_2014.csv", "r")
dict_RNA = dict()
for line in fh:
    gene_name = line.split("|")[0][1:].strip()
    line = fh.readline()
    rna_seq = line.strip()
    dict_RNA[gene_name] = rna_seq
fh.close()


finder = GeneFinder(list(dict_RNA.keys()), list(dict_RNA.values()), species = "mouse", release = 111)
finder._Finder()