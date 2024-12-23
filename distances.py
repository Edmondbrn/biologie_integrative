import os
import pandas as pd
import numpy as np
import ast
import pyensembl as pb
from numba import njit

os.chdir(os.path.dirname(__file__))

class Distances():

    """
    This class is for computing distances between splicing sites and protein fixation site on mRNA
    """


    def __init__(this):
        this.__data_prot : pd.DataFrame = pd.DataFrame()
        this.__data_splicing : dict[ str : pd.DataFrame ] = dict()
    

    def __LoadDataProt(this, file_name : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data of the protein
        file_name : name of the file of the protein
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        this.__data_prot = pd.read_csv(file_name, sep = sep)
        this.__data_prot["start_genomic"] = this.__data_prot["start_genomic"].apply(ast.literal_eval)
        this.__data_prot["end_genomic"] = this.__data_prot["end_genomic"].apply(ast.literal_eval)
        this.__FilterDataProt()

    def __FilterDataProt(this):
        # enlever les lignes avec des str sur la colonne start_ensembl
        this.__data_prot = this.__data_prot[(this.__data_prot["start_ensembl"] != "Not found") | (this.__data_prot["start_ensembl"] != "unknown")]

    
    def __LoadDataSplicing(this, path : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data from alternative splicing
        path : path to the folder containing the files
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        for file in os.listdir(path):
            this.__data_splicing[file.split(".")[0]] = pd.read_csv(f"{path}/{file}", sep = sep)
    
    def __IsDataFrameNull(this, df : pd.DataFrame) -> bool:
        """
        Method to check if a dataframe is null
        """
        return df.shape[0] == 0


    def __CreateDistanceFile(this, df : pd.DataFrame, splice_type : str):
        """
        Method to create a csv file containing the distances
        """
        if os.path.exists("distances") == False:
            os.mkdir("distances")
        df.to_csv(f"distances/distances_{splice_type}.csv", index = False, sep = "\t")

    @staticmethod
    @njit(fastmath = True)
    def __compute_distances(start_genomic_first, end_genomic_last, short_splice, share_splice, long_splice):
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif. 
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour A5SS et A3SS
        """
        distances = np.zeros(8, dtype=np.int64)  # 8 distances
        distances[0] = start_genomic_first - short_splice
        distances[1] = start_genomic_first - share_splice
        distances[2] = end_genomic_last - short_splice
        distances[3] = end_genomic_last - share_splice
        distances[4] = start_genomic_first - long_splice
        distances[5] = start_genomic_first - share_splice
        distances[6] = end_genomic_last - long_splice
        distances[7] = end_genomic_last - share_splice
        return distances
    
    @staticmethod
    @njit(fastmath = True)
    def __compute_distances_RI(start_genomic_first, end_genomic_last, RiStart, RiEnd):
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif. 
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour la rétention d'intron
        """
        distances = np.zeros(4, dtype=np.int64)  # 4 distances

        distances[0] = start_genomic_first - RiStart
        distances[1] = start_genomic_first - RiEnd
        distances[2] = end_genomic_last - RiStart
        distances[3] = end_genomic_last - RiEnd
        return distances
    
    @staticmethod
    @njit(fastmath = True)
    def __compute_distances_SE(start_genomic_first, end_genomic_last, upstreamEnd, DownstreamStart):
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif.
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour l'exon sauté (skipping exon)
        """
        distances = np.zeros(4, dtype=np.int64)
        distances[0] = start_genomic_first - upstreamEnd
        distances[1] = start_genomic_first - DownstreamStart
        distances[2] = end_genomic_last - upstreamEnd
        distances[3] = end_genomic_last - DownstreamStart
        return distances


    def distanceA5SS(this, splice_type: str = "A5SS_+"):
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]]
        results = []  # On stockera ici toutes les lignes calculées
        for i in range(this.__data_prot.shape[0]):
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue

            # Récupération des valeurs start/end (première et dernière coordonnée)
            start_first = this.__data_prot.loc[i, "start_genomic"][0]
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]
            
            for j in range(splicing_same_gene.shape[0]):
                short_sp = splicing_same_gene.iloc[j]["shortSplice"]
                long_sp = splicing_same_gene.iloc[j]["longSplice"]
                share_sp = splicing_same_gene.iloc[j]["shareSplice"]
                
                dist_array = Distances.__compute_distances(start_first, end_last, short_sp, share_sp, long_sp)
                
                # Construire un dict pour tout remettre dans un DataFrame final
                row = {
                    "prot_start_short_splice_start": dist_array[0],
                    "prot_start_short_splice_end":   dist_array[1],
                    "prot_end_short_splice_start":   dist_array[2],
                    "prot_end_short_splice_end":     dist_array[3],
                    "prot_start_long_splice_start":  dist_array[4],
                    "prot_start_long_splice_end":    dist_array[5],
                    "prot_end_long_splice_start":    dist_array[7],
                    "prot_end_long_splice_end":      dist_array[6],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"],
                    "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                results.append(row)
        # Convertir la liste de dicts en DataFrame
        df = pd.DataFrame(results)
        this.__CreateDistanceFile(df, splice_type)  # Suppose qu’on modifie __CreateDistanceFile pour prendre un DataFrame complet

                 
    def distanceA3SS(this, splice_type : str = "A3SS_+") :
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]]
        results = []  # On stockera ici toutes les lignes calculées
        for i in range(this.__data_prot.shape[0]):
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            # Récupération des valeurs start/end (première et dernière coordonnée)
            start_first = this.__data_prot.loc[i, "start_genomic"][0]
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]
            
            for j in range(splicing_same_gene.shape[0]):
                short_sp = splicing_same_gene.iloc[j]["shortSplice"]
                long_sp = splicing_same_gene.iloc[j]["longSplice"]
                share_sp = splicing_same_gene.iloc[j]["shareSplice"]
                dist_array = Distances.__compute_distances(start_first, end_last, short_sp, share_sp, long_sp)
                # Construire un dict pour tout remettre dans un DataFrame final
                row = {
                    "prot_start_short_splice_start": dist_array[1],
                    "prot_start_short_splice_end":   dist_array[0],
                    "prot_end_short_splice_start":   dist_array[3],
                    "prot_end_short_splice_end":     dist_array[2],
                    "prot_start_long_splice_start":  dist_array[5],
                    "prot_start_long_splice_end":    dist_array[4],
                    "prot_end_long_splice_start":    dist_array[6],
                    "prot_end_long_splice_end":      dist_array[7],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"],
                    "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                results.append(row)
        # Convertir la liste de dicts en DataFrame
        df = pd.DataFrame(results)
        this.__CreateDistanceFile(df, splice_type)  # Suppose qu’on modifie __CreateDistanceFile pour prendre un DataFrame complet

    def distanceRI(this, splice_type: str = "RI_+") -> None:
        """
        Fonction pour obtenir les coordonnées des sites de splicing et les distances avec les sites de fixation des protéines
        splice_type : type de splicing

        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "RiStart", "RiEnd"]]  # récupère les coordonnées des sites de splicing
        results = []
        # Browse all the protein fixation sites
        for i in range(this.__data_prot.shape[0]):
            # Select the splicing sites of the same gene to avoid heavy computational loss
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            start_first = this.__data_prot.loc[i, "start_genomic"][0]
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]
            for j in range(splicing_same_gene.shape[0]):
                RiStart = splicing_same_gene.iloc[j]["RiStart"]
                RiEnd = splicing_same_gene.iloc[j]["RiEnd"]
                dist_array = Distances.__compute_distances_RI(start_first, end_last, RiStart, RiEnd)
                row = {
                    "prot_start_RiStart": dist_array[0],
                    "prot_start_RiEnd": dist_array[1],
                    "prot_end_RiStart": dist_array[2],
                    "prot_end_RiEnd": dist_array[3],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"],
                    "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                results.append(row)
        df = pd.DataFrame(results)
        this.__CreateDistanceFile(df, splice_type)
        
    def distanceSE(this, splice_type: str = "SE") -> None:
        """
        Fonction pour obtenir les coordonnées des sites de splicing et les distances avec les sites de fixation des protéines
        splice_type : type de splicing
        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "upstreamEnd", "DownstreamStart"]]
        results = []
        # Browse all the protein fixation sites
        for i in range(this.__data_prot.shape[0]):
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            start_first = this.__data_prot.loc[i, "start_genomic"][0] # extrait les coordonnées de fixation de la protéine
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]
            for j in range(splicing_same_gene.shape[0]):
                upstreamEnd = splicing_same_gene.iloc[j]["upstreamEnd"]
                DownstreamStart = splicing_same_gene.iloc[j]["DownstreamStart"]
                dist_array = Distances.__compute_distances_SE(start_first, end_last, upstreamEnd, DownstreamStart)
                row = {
                    "prot_start_upstreamEnd": dist_array[0],
                    "prot_start_DownstreamStart": dist_array[1],
                    "prot_end_upstreamEnd": dist_array[2],
                    "prot_end_DownstreamStart": dist_array[3],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"],
                    "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                results.append(row)
        df = pd.DataFrame(results)
        this.__CreateDistanceFile(df, splice_type)


    @staticmethod
    @njit(fastmath = True)
    def __compute_distances_MXE(start_genomic_first, end_genomic_last, FirstExonStart, FirstExonEnd, SecondExonStart, SecondExonEnd, upstreamEE, downstreamES):
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif.
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour mutually exclusive exon
        """
        distances = np.zeros(12, dtype=np.int64)  # 16 distances
        distances[0] = start_genomic_first - upstreamEE
        distances[1] = start_genomic_first - FirstExonStart
        distances[2] = start_genomic_first - FirstExonEnd
        distances[3] = start_genomic_first - downstreamES
        distances[4] = end_genomic_last - upstreamEE
        distances[5] = end_genomic_last - FirstExonStart
        distances[6] = end_genomic_last - FirstExonEnd
        distances[7] = end_genomic_last - downstreamES
        distances[8] = start_genomic_first - SecondExonStart
        distances[9] = start_genomic_first - SecondExonEnd
        distances[10] = end_genomic_last - SecondExonStart
        distances[11] = end_genomic_last - SecondExonEnd
        return distances


    def distanceMXE(this, splice_type: str = "MSE") -> None:
            """
            Fonction pour obtenir les coordonnées des sites de splicing et les distances avec les sites de fixation des protéines
            splice_type : type de splicing
            """
            splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "1stExonStart", "1stExonEnd", "2ndExonStart", "2ndExonEnd", "upstreamEE", "downstreamES"]]
            results = []
            # Browse all the protein fixation sites
            for i in range(this.__data_prot.shape[0]):
                splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
                if this.__IsDataFrameNull(splicing_same_gene):
                    continue
                start_first = this.__data_prot.loc[i, "start_genomic"][0] # extrait les coordonnées de fixation de la protéine
                end_last = this.__data_prot.loc[i, "end_genomic"][-1]
                for j in range(splicing_same_gene.shape[0]):
                    FirstExonStart = splicing_same_gene.iloc[j]["1stExonStart"]
                    FirstExonEnd = splicing_same_gene.iloc[j]["1stExonEnd"]
                    SecondExonStart = splicing_same_gene.iloc[j]["2ndExonStart"]
                    SecondExonEnd = splicing_same_gene.iloc[j]["2ndExonEnd"]
                    upstreamEE = splicing_same_gene.iloc[j]["upstreamEE"]
                    downstreamES = splicing_same_gene.iloc[j]["downstreamES"]

                    dist_array = Distances.__compute_distances_MXE(start_first, end_last, FirstExonStart, FirstExonEnd, SecondExonStart, SecondExonEnd, upstreamEE, downstreamES)
                    row = {
                        "prot_start_upstreamEE": dist_array[0],
                        "prot_start_FirstExonStart": dist_array[1],
                        "prot_start_FirstExonEnd": dist_array[2],
                        "prot_start_downstreamES": dist_array[3],
                        "prot_end_upstreamEE": dist_array[4],
                        "prot_end_FirstExonStart": dist_array[5],
                        "prot_end_FirstExonEnd": dist_array[6],
                        "prot_end_downstreamES": dist_array[7],
                        "prot_start_SecondExonStart": dist_array[8],
                        "prot_start_SecondExonEnd": dist_array[9],
                        "prot_end_SecondExonStart": dist_array[10],
                        "prot_end_SecondExonEnd": dist_array[11],
                        "transcript_ID": this.__data_prot.loc[i, "ensembl_id"],
                        "prot_seq": this.__data_prot.loc[i, "seq"]
                    }
                    results.append(row)
            df = pd.DataFrame(results)
            this.__CreateDistanceFile(df, splice_type)

    def ComputeDistance(this, splice_type : str = ""):
        
        distance_calculator = {
            "A5SS_+": this.distanceA5SS,
            "A5SS_-": this.distanceA3SS,
            "A3SS_+": this.distanceA3SS,
            "A3SS_-": this.distanceA5SS,
            "RI_-": this.distanceRI ,
            "RI_+": this.distanceRI ,
            "SE_-": this.distanceSE ,
            "SE_+": this.distanceSE ,
            "MXE_-": this.distanceMXE,
            "MXE_+": this.distanceMXE,
            }

        return distance_calculator[splice_type](splice_type)
        
    def start(this, path_prot : str, path_splicing : str, splice_type : str):
        """
        Method to start the computation of the distances
        path_prot : path to the file containing the protein data with the csv file at the end
        path_splicing : path to the folder containing the splicing data
        """
        this.__LoadDataProt(path_prot)
        this.__LoadDataSplicing(path_splicing)
        this.ComputeDistance(splice_type)

if __name__ == "__main__":
    dist = Distances()
    splice_types = ["A5SS_+", "A5SS_-", "A3SS_+", "A3SS_-", "RI_+", "RI_-", "SE_-", "SE_+", "MXE_-", "MXE_+"]
    for splice_type in splice_types:
        dist.start(path_prot = "data_filteredFMRP.tsv", path_splicing = "filteredRmats", splice_type = splice_type)
        print(f"Distances for {splice_type} computed")