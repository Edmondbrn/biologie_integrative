import os
import pytest
import sys
import pandas as pd
from unittest.mock import patch, MagicMock
from PyQt6.QtCore import QThread
from PyQt6.QtTest import QSignalSpy

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Remplacez mon_package.mon_module par le chemin réel de votre code
from src.Scripts.Back.DistanceWorkerAll import ParallelDistancesWorkerAll

@pytest.fixture
def df_ref():
    """
    DataFrame de référence minimal pour df_ref.
    """
    data = {
        "ensembl_id": ["ENST000001", "ENST000002"],
        "GeneID": ["GENE1", "GENE2"],
        "seq": ["MASEQ1", "MASEQ2"],
        "start_genomic": [100, 200],
    }
    return pd.DataFrame(data)

@pytest.fixture
def input_dfs():
    """
    Dictionnaire de DataFrames pour différents splices.
    """
    data_splice1 = {
        "GeneID": ["GENE1", "GENE2", "GENE1"],
        "coord_B": [155, 255, 160],
    }
    data_splice2 = {
        "GeneID": ["GENE1", "GENE2"],
        "coord_B": [300, 400],
    }
    return {
        "splice1": pd.DataFrame(data_splice1),
        "splice2": pd.DataFrame(data_splice2),
    }

@pytest.fixture
def comparison_couples():
    """
    Par exemple, un dict de listes de tuples pour chaque 'splice'.
    """
    return {
        "splice1": [("start_genomic", "coord_B")],
        "splice2": [("start_genomic", "coord_B")],
    }

@pytest.fixture
def mock_bdd():
    """
    Mock de l'objet EnsemblRelease (pyensembl).
    """
    transcript_mock = MagicMock()
    transcript_mock.exon_intervals = [(100, 120), (130, 150)]
    bdd_mock = MagicMock()
    bdd_mock.transcript_by_id.return_value = transcript_mock
    return bdd_mock

@pytest.fixture
def output_dir(tmp_path):
    """
    Utilisation du dossier temporaire de pytest.
    """
    return str(tmp_path / "output")

@pytest.fixture
def mock_parallel_start_manual_all():
    """
    Patch de la fonction parallel_start_manual_all pour éviter l'exécution réelle
    et capturer ses appels.
    """
    # IMPORTANT : Remplacez mon_package.mon_module par le chemin exact
    # où parallel_start_manual_all est importé dans ParallelDistancesWorkerAll.
    with patch("src.Scripts.Back.distances_utils.parallel_start_manual_all") as mock_func:
        yield mock_func

def test_parallel_distances_worker_all(
    df_ref,
    input_dfs,
    comparison_couples,
    mock_bdd,
    output_dir,
    mock_parallel_start_manual_all
):
    """
    Test unitaire de la classe ParallelDistancesWorkerAll.
    On vérifie :
    - Que parallel_start_manual_all est bien appelée.
    - Que les signaux progress_changed et finished_signal sont émis.
    """

    # Instanciation de la classe
    worker = ParallelDistancesWorkerAll(
        df_ref=df_ref,
        input_dfs=input_dfs,
        comparison_couples=comparison_couples,
        bdd=mock_bdd,
        output_dir=output_dir,
        file_basename="test_parallel",
        n_processes=2
    )

    # On utilise QSignalSpy pour capturer les signaux
    progress_spy = QSignalSpy(worker.progress_changed)
    finished_spy = QSignalSpy(worker.finished_signal)

    # Lancement synchrone de run() (au lieu de .start()) pour simplifier le test
    worker.run()

    # Vérification : parallel_start_manual_all a bien été appelée 1 fois
    assert mock_parallel_start_manual_all.call_count == 1, (
        "La fonction parallel_start_manual_all n'a pas été appelée exactement une fois." 
    )

    # On peut vérifier les arguments passés à parallel_start_manual_all
    called_args, called_kwargs = mock_parallel_start_manual_all.call_args
    assert "df_ref" in called_kwargs, "df_ref non présent dans les arguments appelés"
    assert "input_dfs" in called_kwargs, "input_dfs non présent dans les arguments appelés"
    assert called_kwargs["df_ref"].equals(df_ref), "df_ref passé en argument est différent de l'original."
    assert called_kwargs["input_dfs"] == input_dfs, "input_dfs passé en argument est différent de l'original."
    assert called_kwargs["output_dir"] == output_dir, "output_dir passé en argument est incorrect."
    assert called_kwargs["output_basename"] == "test_parallel", "file_basename passé en argument est incorrect."
    assert called_kwargs["n_cores"] == 2, "n_cores passé en argument est incorrect."

    # Progression : vu qu'on a mocké parallel_start_manual_all, la callback interne
    # ne fera peut-être rien. Mais on peut simuler un renvoi "manuel" si nécessaire.
    # Pour l'instant, on vérifie juste que le signal n'a pas été émis
    # (à moins que votre implémentation d'origine déclenche la callback).
    assert len(progress_spy) == 0, "Des progress_changed ont été émis alors qu'il n'y a pas eu de callback."

    # Vérification du signal finished_signal : la classe l'émet à la fin
    assert len(finished_spy) == 1, "Le signal finished_signal n'a pas été émis exactement une fois."
