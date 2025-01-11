from PyQt6.QtCore import QThread, pyqtSignal
import pandas as pd
import pyensembl as pb
import numpy as np
import os
from .distances_utils import ComputeDistanceManuel_wrapper, parallel_start_manual, FilterDataProt, fill_rna_row


class DistancesWorker(QThread):
    # On déclare les signaux dans la classe
    progress_changed = pyqtSignal(int)
    finished_signal = pyqtSignal()
    error_signal = pyqtSignal(str)

    def __init__(self,
                df_ref,
                df_second,
                comparison_couples,
                output_dir,
                release: int,
                species: str,
                file_basename="distances"):
        super().__init__()
        self.df_ref = df_ref
        self.df_second = df_second
        self.comparison_couples = comparison_couples
        self.output_dir = output_dir
        # On NE crée PAS l'objet bdd ici, on stocke seulement les infos
        self.release = release
        self.species = species
        self.file_basename = file_basename

    def run(self):
        """
        Méthode principale lancée quand on fait .start() sur le QThread.
        -> C'est ici, dans le "vrai" thread worker,
           qu'on doit instancier l'objet EnsemblRelease.
        """
        # Création (et téléchargement/ index) dans le thread *courant*
        self.bdd = pb.EnsemblRelease(release=self.release, species=self.species)
        self.bdd.download()
        self.bdd.index()

        # Maintenant, l'accès à self.bdd se fait dans le même thread => OK
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
            try: # si l'id n'est pas connu dans la base de données
                transcript : pb.Transcript = self.bdd.transcript_by_id(row_ref["ensembl_id"])
            except Exception as e:
                continue
            exon_pos_list = transcript.exon_intervals
            # filtre les données pour limiter les calculs
            df_same_gene : pd.DataFrame = data_splicing.loc[data_splicing["GeneID"] == row_ref["GeneID"]]
            for y in range(len(df_same_gene)): # parcours des différentes lignes du 2e tableau
                row_compare = df_same_gene.iloc[y]
                idx_couple = []
                for couple in self.comparison_couples: # on parcours les différentes combinaisons à effectuer
                    try:
                        array_coord = np.array([int(row_ref[couple[0]]), int(row_compare[couple[1]])])
                        idx_couple.append(array_coord)
                    except Exception as e:
                        self.error_signal.emit(f"Error while converting coordinates into integers. Please check your data and reload the window: {e}")
                        return
                idx_couple = np.array(idx_couple) # conversion en matrice numpy pour les performances
                # Calcul des distances ADN et ARN
                dist_array, flag_array, err_message_array = ComputeDistanceManuel_wrapper(idx_couple, exon_pos_list)
                
                self.progress_changed.emit(i+1)  # émettre le signal de progression
                row_dna = {"transcript_ID": row_ref.get("ensembl_id", ""), "prot_seq": row_ref.get("seq","")}
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
                    row_ref.get("ensembl_id", ""),
                    row_ref.get("seq", "")
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

class ParallelDistancesWorker(QThread):
    progress_changed = pyqtSignal(int)
    finished_signal = pyqtSignal()

    def __init__(self,
                 df_ref : pd.DataFrame,
                 df_splicing : pd.DataFrame,
                 comparison_couples : list[tuple[str, str]],
                 release : int,
                 species : str,
                 output_dir : str,
                 file_basename : str,
                 n_processes : int = 4,
                 parent=None):
        
        super().__init__(parent)
        self.df_ref = df_ref
        self.df_splicing = df_splicing
        self.comparison_couples = comparison_couples
        self.processes = n_processes
        self.release = release
        self.species = species
        self.output_dir = output_dir
        self.file_basename = file_basename

    def run(self):
        self.bdd = pb.EnsemblRelease(release = self.release, species = self.species)
        self.bdd.download()
        self.bdd.index()
        # Définir la fonction callback que l’on passera à parallel_start_manual
        def progress_callback(rows_done: int):
            # Emettre le signal PyQt
            self.progress_changed.emit(rows_done)

        # Lancer la fonction
        parallel_start_manual(
            df_ref=self.df_ref,
            df_splicing=self.df_splicing,
            comparison_couples=self.comparison_couples,
            bdd=self.bdd,
            output_dir=self.output_dir,
            output_basename=self.file_basename,
            n_cores= self.processes ,  # ou ce que vous voulez
            progress_callback=progress_callback
        )

        # Une fois fini :
        self.finished_signal.emit()

if __name__ == "__main__":
    df_ref = pd.read_csv(".\src\\Ressources\\data\\data_filteredfinal.tsv", sep = "\t")
    df_second = pd.read_csv(".\src\\Ressources\\filteredRmats\A5SS_+.csv")
    comparison_couples = [("start_genomic", "shortSplice"), ("end_genomic", "longSplice")]
    output_dir = ".\src\\Ressources\\output_calcul"
    release = 102
    species = "mus_musculus"
    dist = DistancesWorker(df_ref= df_ref,
                           df_second= df_second,
                           comparison_couples= comparison_couples,
                           output_dir= output_dir,
                           release= release,
                           species= species)
    dist.start()