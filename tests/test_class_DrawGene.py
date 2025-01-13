import sys
import os
# ajout de la racine dans le chemin pour pouvoir importer les modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from src.Scripts.Back.DrawGene import GeneImage

@pytest.fixture
def example_exon_pos():
    return [
        [97543299, 97544702], [97547885, 97548026], [97564044, 97564188],
        [97658624, 97658804], [97700407, 97700550], [97770814, 97770934]
    ]

@pytest.fixture
def example_marker_pos():
    return [97647885, 97770934]

@pytest.fixture
def example_exon_pos_false():
    return [
        ["qd", 97544702], ["szef", 97548026], [97564044, 97564188],
        [97658624, 97658804], [97700407, "sedf"], [97770814, 97770934]
    ]



def test_geneimage_draw(example_exon_pos, example_marker_pos):
    """
    Test si la méthode draw de la classe GeneImage fonctionne correctement (vérification visuelle du graphique)
    """
    gene = GeneImage(
        exon_intervals=example_exon_pos,
        marker_pos=example_marker_pos
    )
    assert gene.show() == None 

def test_geneimage_save(example_exon_pos, example_marker_pos):
    """
    Vérification de la sauvegarde de l'image
    """
    gene = GeneImage(
        exon_intervals=example_exon_pos,
        marker_pos=example_marker_pos
    )
    assert gene.save("test.png") == None 

def test_geneimage_border_arguments(example_exon_pos, example_marker_pos):
    """
    Vérification de l'ajout de bordures
    """
    gene = GeneImage(
        exon_intervals=example_exon_pos,
        marker_pos=example_marker_pos,
        bar_xmin=97647885,
        bar_xmax=97648785
    )
    # Vérification des positions des bordures
    assert gene.barColorXmin < gene.barColorXmax
    # vérifcation visuelle
    assert gene.show() == None

    gene2 = GeneImage(
        exon_intervals=example_exon_pos,
        marker_pos=example_marker_pos,
        bar_xmin=97698785,
        bar_xmax=97647885
    )
    # test si border min est plus grand que border max
    assert gene2.barColorXmin > gene2.barColorXmax
    # vérifcation visuelle
    assert gene2.show() == None

def test_geneimage_incorrect_border_arguments(example_exon_pos, example_marker_pos):
    """
    Vérification de la gestion incorrecte des arguments de bordure (envoie une valueError pour la GUI)
    """
    with pytest.raises(ValueError):
        gene = GeneImage(
            exon_intervals=example_exon_pos,
            marker_pos=example_marker_pos,
            bar_xmin="ézez",
            bar_xmax=97647885
        )
        gene.show()

    with pytest.raises(ValueError):
        gene = GeneImage(
            exon_intervals=example_exon_pos,
            marker_pos=example_marker_pos,
            bar_xmin=97647885,
            bar_xmax="ézez"
        )
        gene.show()
    
def test_geneimage_incorrect_exon_intervals_arguments(example_exon_pos_false, example_marker_pos):
    """
    Vérification de la gestion incorrecte des arguments de exon_intervals (envoie une valueError pour la GUI)
    """
    with pytest.raises(ValueError):
        gene = GeneImage(
            exon_intervals=example_exon_pos_false,
            marker_pos=example_marker_pos
        )
        gene.show()

    with pytest.raises(ValueError):
        gene = GeneImage(
            exon_intervals=example_exon_pos_false,
            marker_pos=example_marker_pos
        )
        gene.exonIntervals = "ez"
        gene.show()

