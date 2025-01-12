import os
import sys
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from PyQt6.QtTest import QSignalSpy


from src.Scripts.Back.DistanceWorkerAll import DistancesWorkerAll

# =========================
# Fixtures de test
# =========================

@pytest.fixture
def df_ref():
    """
    DataFrame de référence minimal.
    """
    data = {
        "ensembl_id": ["ENST000001", "ENST000002"],
        "GeneID": ["GENE1", "GENE2"],
        "seq": ["MASEQ1", "MASEQ2"],
        # Exemples de colonnes supplémentaires (coord_A, coord_B, etc.)
        "start_genomic": [100, 200],
        "end_genomic": [150, 250],
    }
    return pd.DataFrame(data)

@pytest.fixture
def input_df():
    """
    Dictionnaire contenant des DataFrames d'exemples pour différentes conditions de 'splice'.
    On crée volontairement des lignes correspondant aux GeneID du df_ref pour les tester.
    """
    data_splice1 = {
        "GeneID": ["GENE1", "GENE2", "GENE1"],
        "coord_B": [155, 255, 160],  # Valeurs fictives
    }
    data_splice2 = {
        "GeneID": ["GENE1", "GENE2"],
        "coord_B": [300, 400],  # Valeurs fictives
    }
    return {
        "splice1": pd.DataFrame(data_splice1),
        "splice2": pd.DataFrame(data_splice2),
    }

@pytest.fixture
def comparison_couples():
    """
    Exemple de couples à comparer.
    """
    # Par exemple, on compare la position "coord_A" du df_ref avec "coord_B" des DataFrames d'entrée
    return {
        "splice1": [("start_genomic", "coord_B")],
        "splice2": [("start_genomic", "coord_B")],
    }

@pytest.fixture
def mock_bdd():
    """
    Mock de l'objet EnsemblRelease (ou bdd) pour éviter de charger une vraie base.
    On mock la méthode transcript_by_id pour qu'elle renvoie un objet transcript minimal.
    """
    transcript_mock = MagicMock()
    # On définit quelques exons fictifs
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
# Mocks des fonctions
# =========================

@pytest.fixture
def mock_filter_data_prot():
    """
    Mock de FilterDataProt qui retourne simplement le DataFrame sans modification,
    ou filtre en fonction de votre logique de test.
    """
    with patch("src.Scripts.Back.distances_utils.FilterDataProt") as mock_func:
        mock_func.side_effect = lambda x: x  # on renvoie x sans le modifier
        yield mock_func

@pytest.fixture
def mock_compute_distance_manual():
    """
    Mock de ComputeDistanceManual qui retourne des distances fictives.
    """
    with patch("src.Scripts.Back.distances_utils.ComputeDistanceManual") as mock_func:
        # Par exemple, on renvoie:
        # dist_array = [10, 20, ...], flag_array = [False, ...], err_message_array = ["", ...]
        mock_func.return_value = (
            np.array([10]),      # dist_array
            np.array([False]),   # flag_array
            np.array([""])       # err_message_array
        )
        yield mock_func

@pytest.fixture
def mock_fill_rna_row():
    """
    Mock de la fonction fill_rna_row qui construit une ligne pour le DataFrame RNA.
    """
    with patch("src.Scripts.Back.distances_utils.fill_rna_row") as mock_func:
        # On renvoie un dict minimal ; vous pouvez ajouter plus de champs si nécessaire
        mock_func.side_effect = lambda rna_indices, dist_array, flag_array, err_message_array, transcript_id, seq: {
            "transcript_ID": transcript_id,
            "prot_seq": seq,
            "some_rna_field": "rna_value"
        }
        yield mock_func


# =========================
# Test principal
# =========================

@pytest.mark.usefixtures(
    "mock_filter_data_prot",
    "mock_compute_distance_manual",
    "mock_fill_rna_row"
)
def test_distances_worker_all(
    df_ref,
    input_df,
    comparison_couples,
    mock_bdd,
    output_dir
):
    """
    Test de la méthode principale 'run' (appelée par 'start_manual_all').
    """


    # Instanciation de la classe
    worker = DistancesWorkerAll(
        df_ref=df_ref,
        input_df=input_df,
        comparison_couples=comparison_couples,
        output_dir=output_dir,
        release = 102,
        species = "mus_musculus",
        file_basename="test_distances"
    )

    worker.bdd = mock_bdd  # On remplace la vraie base par le mock

    # On utilise QSignalSpy pour observer les signaux
    progress_spy = QSignalSpy(worker.progress_changed)
    finished_spy = QSignalSpy(worker.finished_signal)

    # Lancement du thread en mode synchrone (pour simplifier les tests)
    # Normalement, on ferait worker.start() et on gérerait l'asynchrone,
    # mais dans un test unitaire, on peut appeler 'run()' directement.
    worker.run()

    # Vérification : est-ce que la progression a été émise au moins une fois ?
    assert len(progress_spy) > 0, "Le signal de progression n'a pas été émis."

    # Vérification : est-ce que le signal de fin a été émis exactement 1 fois ?
    assert len(finished_spy) == 1, "Le signal de fin n'a pas été émis une seule fois."
    # Vérification de la création des fichiers de sortie
    # Deux conditions 'splice1' et 'splice2' => 2 paires de fichiers (dna_..., rna_...)
    expected_files = [
        f"dna_splice1_test_distances.csv",
        f"rna_splice1_test_distances.csv",
        f"dna_splice2_test_distances.csv",
        f"rna_splice2_test_distances.csv",
    ]
    for fname in expected_files:
        outpath = os.path.join(output_dir, fname)
        print(outpath)
        assert os.path.exists(outpath), f"Le fichier de sortie {fname} n'a pas été créé."

    # On charge les fichiers CSV et vérifier leur contenu
    # pour s'assurer qu'on y trouve les colonnes attendues ou la bonne taille.:
    dna_splice1 = pd.read_csv(os.path.join(output_dir, "dna_splice1_test_distances.csv"), sep="\t")
    assert "transcript_ID" in dna_splice1.columns, "Colonne 'transcript_ID' manquante dans dna_splice1"
