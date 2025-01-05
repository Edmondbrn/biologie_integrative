from PyQt6.QtCore import QThread, pyqtSignal
from distances_utils import FilterDataProt, fill_rna_row
import pandas as pd
import pyensembl as pb
import numpy as np
import os
from distances_utils import ComputeDistanceManual

class DistancesWorker(QThread):
    # On déclare les signaux dans la classe
    progress_changed = pyqtSignal(int)
    finished_signal = pyqtSignal()

    def __init__(self, df_ref, df_second, comparison_couples, output_dir, bdd : pb.EnsemblRelease, file_basename="distances"):
        super().__init__()
        self.df_ref = df_ref
        self.df_second = df_second
        self.comparison_couples = comparison_couples
        self.output_dir = output_dir
        self.bdd = bdd
        self.file_basename = file_basename

    def run(self):
        """
        Méthode principale qui sera lancée quand on fait .start() sur le thread.
        """
        # Placez votre logique "lourde" ici ou appelez une fonction
        self.start_manual()

    def start_manual(self) -> None:
        """
        Méthode principale pour lancer tous les calculs de distances 'manuelles'.
        """
        data_prot : pd.DataFrame = self.df_ref
        data_prot = FilterDataProt(data_prot) # enlève les lignes avec des valeurs manquantes
        data_splicing : pd.DataFrame = self.df_second
        results_dna, results_rna = [], []
        nb_couple = len(self.comparison_couples)
        for i in range(len(data_prot)): # parcours des différentes lignes de la table de fixation des protéines
            row_ref = data_prot.iloc[i] # stocke la ligne pour la lisibilité
            # on récupère les coordonnées des exons
            transcript : pb.Transcript = self.bdd.transcript_by_id(row_ref["ensembl_id"])
            exon_pos_list = transcript.exon_intervals
            # filtre les données pour limiter les calculs
            df_same_gene : pd.DataFrame = data_splicing.loc[data_splicing["GeneID"] == row_ref["GeneID"]]
            for y in range(len(df_same_gene)): # parcours des différentes lignes du 2e tableau
                row_compare = df_same_gene.iloc[y]
                idx_couple = []
                for couple in self.comparison_couples: # on parcours les différentes combinaisons à effectuer
                    array_coord = np.array([int(row_ref[couple[0]]), int(row_compare[couple[1]])])
                    idx_couple.append(array_coord)
                idx_couple = np.array(idx_couple) # conversion en matrice numpy pour les performances
                # Calcul des distances ADN et ARN
                dist_array, flag_array, err_message_array = ComputeDistanceManual(idx_couple, exon_pos_list)
                
                self.progress_changed.emit(i+1)  # émettre le signal de progression
                
                row_dna = {"transcript_ID": row_ref["ensembl_id"], "prot_seq": row_ref["seq"]}
                rna_indices = {}
                for y, couple in enumerate(self.comparison_couples):
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
        df_dna.to_csv(f"{self.output_dir}/dna_{self.file_basename}.csv", sep="\t", index=False)
        df_rna.to_csv(f"{self.output_dir}/rna_{self.file_basename}.csv", sep="\t", index=False)
        self.finished_signal.emit()  # émettre le signal de fin

