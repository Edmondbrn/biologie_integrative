from PyQt6.QtCore import QThread, pyqtSignal
import pandas as pd
import pyensembl as pb
import numpy as np
import os
from .distances_utils import ComputeDistanceManual, parallel_start_manual_all, FilterDataProt, fill_rna_row


class DistancesWorkerAll(QThread):
    # On déclare les signaux dans la classe
    progress_changed = pyqtSignal(int)
    finished_signal = pyqtSignal()

    def __init__(self, df_ref : pd.DataFrame, 
                 input_df : dict[str : pd.DataFrame], 
                 comparison_couples : dict[str : [tuple[str, str]]], 
                 output_dir : str, 
                 bdd : pb.EnsemblRelease, 
                 file_basename="distances"):
        
        super().__init__()
        self.df_ref = df_ref
        self.input_df = input_df
        self.dict_comparison_couples = comparison_couples
        self.output_dir = output_dir
        self.bdd = bdd
        self.file_basename = file_basename

    def run(self):
        """
        Méthode principale qui sera lancée quand on fait .start() sur le thread.
        """
        self.start_manual_all()

    def start_manual_all(self) -> None:
        """
        Méthode principale pour lancer tous les calculs de distances 'manuelles'.
        """
        data_prot : pd.DataFrame = self.df_ref
        data_prot = FilterDataProt(data_prot) # enlève les lignes avec des valeurs manquantes
        cpt = 0
        for splice, data_splicing in self.input_df.items():
            results_dna, results_rna = [], []
            comparison_couples = self.dict_comparison_couples[splice]
            nb_couple = len(comparison_couples)
            for i in range(len(data_prot)): # parcours des différentes lignes de la table de fixation des protéines
                cpt += 1
                row_ref = data_prot.iloc[i] # stocke la ligne pour la lisibilité
                # on récupère les coordonnées des exons
                transcript : pb.Transcript = self.bdd.transcript_by_id(row_ref["ensembl_id"])
                exon_pos_list = transcript.exon_intervals
                # filtre les données pour limiter les calculs
                df_same_gene : pd.DataFrame = data_splicing.loc[data_splicing["GeneID"] == row_ref["GeneID"]]
                for y in range(len(df_same_gene)): # parcours des différentes lignes du 2e tableau
                    row_compare = df_same_gene.iloc[y]
                    idx_couple = []
                    for couple in comparison_couples: # on parcours les différentes combinaisons à effectuer
                        array_coord = np.array([int(row_ref[couple[0]]), int(row_compare[couple[1]])])
                        idx_couple.append(array_coord)
                    idx_couple = np.array(idx_couple) # conversion en matrice numpy pour les performances
                    # Calcul des distances ADN et ARN
                    dist_array, flag_array, err_message_array = ComputeDistanceManual(idx_couple, exon_pos_list)
                    self.progress_changed.emit(cpt)  # émettre le signal de progression
                    row_dna = {"transcript_ID": row_ref["ensembl_id"], "prot_seq": row_ref["seq"]}
                    rna_indices = {}
                    for y, couple in enumerate(comparison_couples):
                        row_dna[f"coord_{couple[0]}"] = row_ref[couple[0]]
                        row_dna[f"coord_{couple[1]}"] = row_compare[couple[1]]
                        row_dna[f"{couple[0]}-{couple[1]}"] = dist_array[y] # construction des lignes pour le DataFrame final de l'ADN
                        rna_indices[nb_couple + y] = f"{couple[0]}-{couple[1]}" # construction des indices pour les distances ARN

                    row_rna = fill_rna_row(
                        rna_indices,
                        dist_array,
                        flag_array,
                        err_message_array,
                        row_ref["ensembl_id"],
                        row_ref["seq"]
                    )
                    results_dna.append(row_dna)
                    results_rna.append(row_rna)
            df_dna = pd.DataFrame(results_dna)
            df_rna = pd.DataFrame(results_rna)
            if not os.path.isdir(self.output_dir):
                os.makedirs(self.output_dir)
            df_dna.to_csv(f"{self.output_dir}/dna_{splice}_{self.file_basename}.csv", sep="\t", index=False)
            df_rna.to_csv(f"{self.output_dir}/rna_{splice}_{self.file_basename}.csv", sep="\t", index=False)
        self.finished_signal.emit()  # émettre le signal de fin


class ParallelDistancesWorkerAll(QThread):
    progress_changed = pyqtSignal(int)
    finished_signal = pyqtSignal()

    def __init__(self,
                 df_ref : pd.DataFrame,
                 input_dfs : dict[str : pd.DataFrame],
                 comparison_couples : dict[str : list[tuple[str, str]]],
                 bdd : pb.EnsemblRelease,
                 output_dir : str,
                 file_basename : str,
                 n_processes : int = 4,
                 parent=None):
        
        super().__init__(parent)
        self.df_ref = df_ref
        self.input_dfs = input_dfs
        self.comparison_couples = comparison_couples
        self.bdd = bdd
        self.processes = n_processes
        self.output_dir = output_dir
        self.file_basename = file_basename

    def run(self):
        # Définir la fonction callback que l’on passera à parallel_start_manual
        def progress_callback(rows_done: int):
            # Emettre le signal PyQt
            self.progress_changed.emit(rows_done)

        # Lancer la fonction
        parallel_start_manual_all(
            df_ref=self.df_ref,
            input_dfs=self.input_dfs,
            comparison_couples=self.comparison_couples,
            bdd=self.bdd,
            output_dir=self.output_dir,
            output_basename=self.file_basename,
            n_cores= self.processes ,  # ou ce que vous voulez
            progress_callback=progress_callback
        )

        # Une fois fini :
        self.finished_signal.emit()