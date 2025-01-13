import sys
import os
import pandas as pd
# ajout de la racine dans le chemin pour pouvoir importer les modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
import src.Scripts.Back.distances_utils as du
import numpy as np


@pytest.fixture
def get_dataframe():
    return pd.read_csv("../src/Ressources/data/data_filteredfinal.tsv", sep="\t")

@pytest.fixture
def get_exon_list():
    return [
        (97543299, 97544702), (97547885, 97548026), (97564044, 97564188),
        (97658624, 97658804), (97700407, 97700550), (97770814, 97770934)
    ]

@pytest.fixture
def get_prot_coord():
    return 97547895

@pytest.fixture
def get_splice_coord_before_prot():
    return 97547885

@pytest.fixture
def get_splice_coord_after_prot():
    return 97548026

@pytest.fixture
def get_splice_coord_next_exon():
    return 97564045

@pytest.fixture
def get_splice_coord_prev_exon():
    return 97544701

@pytest.fixture
def get_splice_coord_next_intron():
    return 97548028


def test_filterdataprot(get_dataframe):
    """
    Test de la fonction filterdataprot (regarde si le nombre de lignes diminue)
    """
    df : pd.DataFrame = get_dataframe
    df_filtered = du.FilterDataProt(df)
    # vérifie si la sortie est un dataframe
    assert isinstance(df_filtered, pd.DataFrame)
    assert df_filtered.shape[0] <= df.shape[0]

    # test si on ne lui donne pas un dataframe en entrée
    with pytest.raises(TypeError):
        df_false = du.FilterDataProt("df")


def test_get_intron_coord(get_exon_list):
    """
    Test de la fonction get_intron_coord (vérifie si les coordonnées des introns sont correctes)
    """
    exon_list = get_exon_list
    intron_list = du.get_intron_coord(exon_list)
    print(intron_list)
    # vérification si la sortie est une liste
    assert isinstance(intron_list, list)
    # vérification si la longueur de la liste est correcte (introns = exons - 1)
    assert len(intron_list) == len(exon_list) - 1

    # test si on ne lui donne pas une liste en entrée
    with pytest.raises(TypeError):
        intron_list = du.get_intron_coord("exon_list")
    
    # test si les exons ont plus de 3 éléments
    with pytest.raises(ValueError):
        intron_list = du.get_intron_coord([(97543299, 97544702, 97544702), (97547885, 97548026,1), (97564044, 97564188,2)])


def test_get_intron_length(get_exon_list):
    """
    Test de la fonction get_intron_length (vérifie si les longueurs des introns sont correctes)
    """
    exon_list = get_exon_list
    intron_list = du.get_intron_coord(exon_list)
    intron_length = du.get_intron_length(intron_list[0][0], intron_list[0][1])
    print(f"Results : {intron_length}")
    print(f"Expected : {intron_list[0][1] - intron_list[0][0]}")
    # vérification si la sortie est un int
    assert isinstance(intron_length, int)
    # vérification si la longueur de l'intron est correcte
    assert intron_length == intron_list[0][1] - intron_list[0][0]

    # test si on ne lui donne pas de int en entrée
    with pytest.raises(TypeError):
        intron_length = du.get_intron_length("intron_list", "5")

def test_convert_dna_to_rna(get_prot_coord, get_splice_coord_after_prot, get_splice_coord_before_prot, get_splice_coord_next_exon, get_splice_coord_next_intron, get_splice_coord_prev_exon, get_exon_list):
    """
    Test des fonctions nécessaires pour convertir les coordonnées de l'ADN en ARN
    """
    prot_coord = get_prot_coord
    splice_coord_after_prot = get_splice_coord_after_prot
    splice_coord_before_prot = get_splice_coord_before_prot
    splice_coord_next_exon = get_splice_coord_next_exon
    splice_coord_next_intron = get_splice_coord_next_intron
    splice_coord_prev_exon = get_splice_coord_prev_exon
    exon_list = get_exon_list


    print(f"Prot coord : {prot_coord}")
    print(f"Splice coord after prot : {splice_coord_after_prot}")
    print(f"Splice coord before prot : {splice_coord_before_prot}")
    print(f"Splice coord next exon : {splice_coord_next_exon}")
    print(f"Splice coord next intron : {splice_coord_next_intron}")
    print(f"Splice coord prev exon : {splice_coord_prev_exon}")

    print(f"Results : {du.convert_dna_to_rna(prot_coord, splice_coord_after_prot, prot_coord - splice_coord_after_prot, exon_list)}")
    assert du.convert_dna_to_rna(prot_coord, splice_coord_after_prot, prot_coord - splice_coord_after_prot, exon_list) == (prot_coord - splice_coord_after_prot, False, 0)

    print(f"Results : {du.convert_dna_to_rna(prot_coord, splice_coord_before_prot, prot_coord - splice_coord_before_prot, exon_list)}")
    assert du.convert_dna_to_rna(prot_coord, splice_coord_before_prot, prot_coord - splice_coord_before_prot, exon_list) == (prot_coord - splice_coord_before_prot, False, 0)

    print(f"Results : {du.convert_dna_to_rna(prot_coord, splice_coord_next_exon, prot_coord - splice_coord_next_exon, exon_list)}")
    assert du.convert_dna_to_rna(prot_coord, splice_coord_next_exon, prot_coord - splice_coord_next_exon, exon_list) == (-132, False, 0)

    print(f"Results : {du.convert_dna_to_rna(prot_coord, splice_coord_next_intron, prot_coord - splice_coord_next_intron, exon_list)}")
    assert du.convert_dna_to_rna(prot_coord, splice_coord_next_intron, prot_coord - splice_coord_next_intron, exon_list) == (prot_coord-splice_coord_next_intron, True, 0)

    print(f"Results : {du.convert_dna_to_rna(prot_coord, splice_coord_prev_exon, prot_coord - splice_coord_prev_exon, exon_list)}")
    assert du.convert_dna_to_rna(prot_coord, splice_coord_prev_exon, prot_coord - splice_coord_prev_exon, exon_list) == (11, False, 0)

    print(f"Resuls : {du.convert_dna_to_rna(prot_coord, 999999999999999, prot_coord - 999999999999, exon_list)}")
    assert du.convert_dna_to_rna(prot_coord, 999999999999999, prot_coord - 999999999999, exon_list) == (0, False, 1)

    print(f"Resuls : {du.convert_dna_to_rna(8888888888888, splice_coord_prev_exon, prot_coord - 8888888888888, exon_list)}")
    assert du.convert_dna_to_rna(8888888888888, splice_coord_prev_exon, prot_coord - 8888888888888, exon_list) == (0, False, 2)

    # tes si prot coord n'est pas un int
    with pytest.raises(TypeError):
        du.convert_dna_to_rna("prot_coord", splice_coord_prev_exon, prot_coord - splice_coord_prev_exon, exon_list)
        du.convert_dna_to_rna(prot_coord, "splice_coord_prev_exon", prot_coord - splice_coord_prev_exon, exon_list)
        du.convert_dna_to_rna(prot_coord, splice_coord_prev_exon, "prot_coord - splice_coord_prev_exon", exon_list)
        du.convert_dna_to_rna(prot_coord, splice_coord_prev_exon, prot_coord - splice_coord_prev_exon, "exon_list")

def test_compute_distances(get_exon_list, get_prot_coord, get_splice_coord_after_prot, get_splice_coord_next_exon):
    """
    Test la fonction compute distance numba (impossible de tester des arguments inccorrects car la fonction est compilée)
    """
    exon_list = get_exon_list
    prot_coord = get_prot_coord
    splice_coord_after_prot = get_splice_coord_after_prot
    splice_coord_next_exon = get_splice_coord_next_exon
    couple = np.array(([prot_coord, splice_coord_after_prot], [prot_coord, splice_coord_next_exon], [splice_coord_next_exon, splice_coord_after_prot]))

    print(f"Results : {du.ComputeDistanceManual(couple, exon_list)}")
    print(f"Exected : {(np.array([  -131, -16150,  16019,   -131,   -132,      1]), np.array([False, False, False, False, False, False]), np.array([0, 0, 0, 0, 0, 0]))}")


def test_warmup_numba():
    """
    Test de la fonction warmup_numba (vérifie si la fonction est bien compilée)
    """
    assert du.warmup_numba() == 0

 


