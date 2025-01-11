import os
import pyensembl as pb

os.chdir(os.path.dirname(__file__))

class Distances():
    """
    This class is for computing distances between splicing sites and protein fixation site on mRNA
    """


    def __init__(this, ensembl_release: int = 102, specy: str = "mus_musculus"):
        this.bdd = pb.EnsemblRelease(release = ensembl_release, species = specy)  # Déclaration de la variable de classe


  