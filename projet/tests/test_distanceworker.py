import os
import pytest
import pandas as pd
import sys
import numpy as np
from unittest.mock import MagicMock, patch

# Si vous utilisez PyQt6
from PyQt6.QtCore import QThread
from PyQt6.QtTest import QSignalSpy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Remplacez mon_package.mon_module par le chemin réel de votre code
from src.Scripts.Back.DistanceWorker import DistancesWorker

@pytest.fixture
def df_ref():
    """
    DataFrame de référence minimal, représentant par exemple des protéines.
    """
    data = {
        "ensembl_id": ["ENST000001", "ENST000002"],
        "GeneID": ["GENE1", "GENE2"],
        "seq": ["MASEQ1", "MASEQ2"],
        "start_genomic": [100, 200],
    }
    return pd.DataFrame(data)

@pytest.fixture
def df_second():
    """
    Second DataFrame, représentant par exemple un tableau de splicing ou autre.
    """
    data = {
        "GeneID": ["GENE1", "GENE2", "GENE1"],
        "coord_B": [155, 255, 160],  # Valeurs fictives
    }
    return pd.DataFrame(data)

@pytest.fixture
def comparison_couples():
    """
    Liste de tuples décrivant les couples de colonnes à comparer.
    """
    # Exemple : On compare coord_A (df_ref) à coord_B (df_second)
    return [("start_genomic", "coord_B")]

@pytest.fixture
def mock_bdd():
    """
    Mock de l'objet EnsemblRelease (pyensembl) pour éviter de charger une vraie base.
    On mock la méthode transcript_by_id pour qu'elle renvoie un objet transcript minimal.
    """
    transcript_mock = MagicMock()
    # On définit quelques intervalles fictifs d'exons
    transcript_mock.exon_intervals = [(100, 120), (130, 150)]

    bdd_mock = MagicMock()
    bdd_mock.transcript_by_id.return_value = transcript_mock
    return bdd_mock

@pytest.fixture
def output_dir(tmp_path):
    """
    Utilise le dossier temporaire fourni par pytest pour les fichiers de sortie.
    """
    return str(tmp_path / "output")

# =========================
# Mocks des fonctions/objets
# =========================

@pytest.fixture
def mock_filter_data_prot():
    """
    Mock de FilterDataProt qui retourne simplement le DataFrame sans modification,
    ou applique une logique simulée si nécessaire.
    """
    with patch("src.Scripts.Back.distances_utils.FilterDataProt") as mock_func:
        mock_func.side_effect = lambda df: df  # On renvoie df tel quel
        yield mock_func

@pytest.fixture
def mock_compute_distance_manual():
    """
    Mock de ComputeDistanceManual qui retourne des valeurs fictives.
    """
    with patch("src.Scripts.Back.distances_utils.ComputeDistanceManual") as mock_func:
        # Par exemple, on renvoie un dist_array = [10], flag_array = [False], err_message_array = [""].
        mock_func.return_value = (
            np.array([10]),
            np.array([False]),
            np.array([""])
        )
        yield mock_func

@pytest.fixture
def mock_fill_rna_row():
    """
    Mock de la fonction fill_rna_row, qui construit une ligne pour le DataFrame RNA.
    """
    with patch("src.Scripts.Back.distances_utils.fill_rna_row") as mock_func:
        # On renvoie un dict minimal ; vous pouvez adapter si nécessaire.
        mock_func.side_effect = lambda rna_indices, dist_array, flag_array, err_message_array, transcript_id, seq: {
            "transcript_ID": transcript_id,
            "prot_seq": seq,
            "some_rna_column": "test_value"
        }
        yield mock_func


# =========================
# Test de la classe DistancesWorker
# =========================

@pytest.mark.usefixtures(
    "mock_filter_data_prot",
    "mock_compute_distance_manual",
    "mock_fill_rna_row"
)
def test_distances_worker(
    df_ref,
    df_second,
    comparison_couples,
    mock_bdd,
    output_dir
):
    """
    Test unitaire de la classe DistancesWorker.
    Vérifie  :
    - Que la progression est signalée via progress_changed
    - Que le signal finished_signal est émis
    - Que les fichiers CSV sont bien créés
    """


    # Instanciation de la classe
    worker = DistancesWorker(
        df_ref=df_ref,
        df_second=df_second,
        comparison_couples=comparison_couples,
        output_dir=output_dir,
        release = 102,
        species = "mus_musculus",
        file_basename="test_distances"
    )
    worker.bdd = mock_bdd  # On remplace l'objet EnsemblRelease par le mock

    # On utilise QSignalSpy pour observer les signaux
    progress_spy = QSignalSpy(worker.progress_changed)
    finished_spy = QSignalSpy(worker.finished_signal)

    # Exécution "synchrone" de run() (au lieu de .start()), pour éviter la gestion asynchrone en test
    worker.run()

    # Vérification : est-ce que la progression a été émise au moins une fois ?
    assert len(progress_spy) > 0, "Le signal de progression n'a pas été émis."

    # Vérification : est-ce que le signal de fin a été émis exactement 1 fois ?
    assert len(finished_spy) == 1, "Le signal de fin n'a pas été émis une seule fois."

    # Vérification de la création des fichiers de sortie (2 fichiers : dna_*.csv et rna_*.csv)
    dna_file = os.path.join(output_dir, "dna_test_distances.csv")
    rna_file = os.path.join(output_dir, "rna_test_distances.csv")

    assert os.path.exists(dna_file), "Le fichier DNA n'a pas été créé."
    assert os.path.exists(rna_file), "Le fichier RNA n'a pas été créé."

    # On lit le contenu des CSV et vérifier par exemple la présence de certaines colonnes
    df_dna = pd.read_csv(dna_file, sep="\t")
    df_rna = pd.read_csv(rna_file, sep="\t")

    assert "transcript_ID" in df_dna.columns, "La colonne 'transcript_ID' est manquante dans le CSV DNA."
    assert "transcript_ID" in df_rna.columns, "La colonne 'transcript_ID' est manquante dans le CSV RNA."