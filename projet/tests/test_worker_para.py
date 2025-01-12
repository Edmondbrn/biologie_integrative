import os
import sys
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from src.Scripts.Back.distances_utils import parallel_start_manual

# On suppose que 'parallel_start_manual' est importé 
# depuis src.Scripts.Back.distances_utils
from src.Scripts.Back.distances_utils import parallel_start_manual

def test_parallel_start_manual_chunking():
    # Création d'un DataFrame df_ref de 10 lignes
    data = {
        "ensembl_id": [f"ENST00000{i}" for i in range(10)],
        "GeneID": [f"Gene{i}" for i in range(10)],
        "seq": [f"Sequence{i}" for i in range(10)]
    }
    df_ref = pd.DataFrame(data)
    
    # DataFrame df_splicing vide ou minimal pour le test
    df_splicing = pd.DataFrame(columns=["GeneID"])
    
    # On se fiche du contenu, juste besoin de passer quelque chose
    comparison_couples = [("start", "end")]
    bdd_mock = MagicMock()  # On va mocker l'objet bdd
    
    # Nombre de cœurs pour chunking prévisible
    n_cores = 3
    # chunk_size = 10 // 3 + 1 = 4
    # => chunk1=4 lignes, chunk2=4 lignes, chunk3=2 lignes

    # Patch la classe 'Pool' importée dans 'distances_utils'
    with patch("src.Scripts.Back.distances_utils.Pool") as MockPool:
        # MockPool est la "classe" mockée.
        # Récupérons l’instance qui sera retournée quand on fera Pool(...)
        mock_pool_instance = MockPool.return_value
        
        # On peut aussi s’assurer que apply_async renvoie immédiatement 
        # un objet simulant la valeur de retour, ex. un MagicMock
        mock_pool_instance.apply_async.return_value = MagicMock()

        def progress_callback_mock(_):
            pass

        def mock_apply_async(func, args=(), callback=None):
            # Simule un résultat
            dummy_df_dna = pd.DataFrame({"dna":[1,2]})
            dummy_df_rna = pd.DataFrame({"rna":[3,4]})
            # Appelle la callback avec ce résultat
            if callback:
                callback((dummy_df_dna, dummy_df_rna))
            return MagicMock()
        mock_pool_instance.apply_async.side_effect = mock_apply_async


        # Appel de la fonction à tester
        parallel_start_manual(
            df_ref=df_ref,
            df_splicing=df_splicing,
            comparison_couples=comparison_couples,
            bdd=bdd_mock,
            output_dir="/tmp",  
            output_basename="test",
            n_cores=n_cores,
            progress_callback=progress_callback_mock
        )

        # Vérifier que apply_async a bien été appelé 3 fois (un par chunk)
        assert mock_pool_instance.apply_async.call_count == 3, (
            f"apply_async devrait être appelé 3 fois (une par chunk), "
            f"mais a été appelé {mock_pool_instance.apply_async.call_count} fois."
        )

        # Récupérer tous les appels (liste de call)
        call_args_list = mock_pool_instance.apply_async.call_args_list

        # Vérifier la taille du chunk passé à chaque appel
        for idx, call_args in enumerate(call_args_list):
            # call_args est un tuple ( (arg1, arg2, ...), {kwarg1: val1, ...} )
            positional_args, keyword_args = call_args
            # positional_args[0] = la fonction => process_chunk
            # positional_args[1]["args"] => arguments passés à process_chunk
            # Mais parfois c’est plus simple de faire :
            args_tuple = keyword_args["args"]  # (chunk, df_splicing, ...)
            chunk_passed = args_tuple[0]       # le 1er argument = chunk

            expected_len = 4 if idx < 2 else 2  # 4, 4, 2
            assert len(chunk_passed) == expected_len, (
                f"Le chunk n°{idx} devrait avoir une taille de {expected_len}, "
                f"mais a {len(chunk_passed)} lignes."
            )
