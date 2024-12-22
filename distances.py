import os
import pandas as pd
import numpy as np
import ast
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


    def __CreateDistanceFile(this, distance : dict, splice_type : str):
        """
        Method to create a csv file containing the distances
        """
        if os.path.exists("distances") == False:
            os.mkdir("distances")
        df = pd.DataFrame(distance)
        df.to_csv(f"distances/distances_{splice_type}.csv", index = False, sep = "\t")

    @njit
    def __compute_distances(start_genomic_first, end_genomic_last, short_splice, share_splice, long_splice):
        """
        Calcule les distances pour un enregistrement. 
        Les arguments sont des entiers qui représentent la première coordonnée 
        d'une liste (start_genomic) et la dernière coordonnée (end_genomic), 
        ainsi que short_splice, share_splice et long_splice.
        """
        distances = np.zeros(8, dtype=np.int64)  # 8 distances
        # prot_start_short_splice_start
        distances[0] = start_genomic_first - short_splice
        # prot_start_short_splice_end
        distances[1] = start_genomic_first - share_splice
        # prot_end_short_splice_start
        distances[2] = end_genomic_last - short_splice
        # prot_end_short_splice_end
        distances[3] = end_genomic_last - share_splice
        # prot_start_long_splice_start
        distances[4] = start_genomic_first - long_splice
        # prot_start_long_splice_end
        distances[5] = start_genomic_first - share_splice
        # prot_end_long_splice_start
        distances[6] = end_genomic_last - long_splice
        # prot_end_long_splice_end
        distances[7] = end_genomic_last - share_splice
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
                
                dist_array = this.__compute_distances(start_first, end_last, short_sp, share_sp, long_sp)
                
                # Construire un dict pour tout remettre dans un DataFrame final
                row = {
                    "prot_start_short_splice_start": dist_array[0],
                    "prot_start_short_splice_end":   dist_array[1],
                    "prot_end_short_splice_start":   dist_array[2],
                    "prot_end_short_splice_end":     dist_array[3],
                    "prot_start_long_splice_start":  dist_array[4],
                    "prot_start_long_splice_end":    dist_array[5],
                    "prot_end_long_splice_start":    dist_array[6],
                    "prot_end_long_splice_end":      dist_array[7],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"],
                    "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                results.append(row)
        
        # Convertir la liste de dicts en DataFrame
        df = pd.DataFrame(results)
        this.__CreateDistanceFile(df, splice_type)  # Suppose qu’on modifie __CreateDistanceFile pour prendre un DataFrame complet

    
    # def distanceA5SS(this, splice_type : str = "A5SS_+"):
    #     splicing_coordinates : pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]] #récupère les coordonnées des sites de splicing
    #     # Browse all the protein fixation sites
    #     for i in range(this.__data_prot.shape[0]):
    #         # Select the splicing sites of the same gene to avoid heavy computationnal loss
    #         splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
    #         if this.__IsDataFrameNull(splicing_same_gene):
    #             continue
    #         # Browse all the splicing sites of the same gene
    #         distance = {"prot_start_short_splice_start" : [], "prot_start_short_splice_end" : [], "prot_end_short_splice_start" : [], "prot_end_short_splice_end" : [],
    #                     "prot_start_long_splice_start" : [], "prot_start_long_splice_end" : [], "prot_end_long_splice_start" : [], "prot_end_long_splice_end" : [],
    #                     "transcript_ID" : [], "prot_seq" : []}
    #         for j in range(splicing_same_gene.shape[0]):
    #             # Compute the distance between the protein fixation site and the splicing site
    #             distance["prot_start_short_splice_start"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["shortSplice"])
    #             distance["prot_start_short_splice_end"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["shareSplice"])
    #             distance["prot_end_short_splice_start"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["shortSplice"])
    #             distance["prot_end_short_splice_end"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["shareSplice"])
    #             distance["prot_start_long_splice_start"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["longSplice"])
    #             distance["prot_start_long_splice_end"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["shareSplice"])
    #             distance["prot_end_long_splice_start"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["longSplice"])
    #             distance["prot_end_long_splice_end"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["shareSplice"])
    #             distance["transcript_ID"].append(this.__data_prot.loc[i, "ensembl_id"])
    #             distance["prot_seq"].append(this.__data_prot.loc[i, "seq"])
    #     this.__CreateDistanceFile(distance, splice_type)
                
                
    def distanceA3SS(this, splice_type : str = "A3SS_+") :
            splicing_coordinates : pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]] #récupère les coordonnées des sites de splicing
            # Browse all the protein fixation sites
            for i in range(this.__data_prot.shape[0]):
                # Select the splicing sites of the same gene to avoid heavy computationnal loss
                splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
                if this.__IsDataFrameNull(splicing_same_gene):
                    continue
                # Browse all the splicing sites of the same gene
                distance = {"prot_start_short_splice_start" : [], "prot_start_short_splice_end" : [], "prot_end_short_splice_start" : [], "prot_end_short_splice_end" : [],
                        "prot_start_long_splice_start" : [], "prot_start_long_splice_end" : [], "prot_end_long_splice_start" : [], "prot_end_long_splice_end" : [],
                        "transcript_ID" : [], "prot_seq" : []}
                
                for j in range(splicing_same_gene.shape[0]):
                # Compute the distance between the protein fixation site and the splicing site
                # "shortSplice", "longSplice", "shareSplice"
                # Compute the distance between the protein fixation site and the splicing site
                    distance["prot_start_short_splice_start"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["shareSplice"])
                    distance["prot_start_short_splice_end"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["shortSplice"])
                    distance["prot_end_short_splice_start"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["shareSplice"])
                    distance["prot_end_short_splice_end"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["shortSplice"])
                    distance["prot_start_long_splice_start"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["shareSplice"])
                    distance["prot_start_long_splice_end"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["longSplice"])
                    distance["prot_end_long_splice_start"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["shareSplice"])
                    distance["prot_end_long_splice_end"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["longSplice"])
                    distance["transcript_ID"].append(this.__data_prot.loc[i, "ensembl_id"])
                    distance["prot_seq"].append(this.__data_prot.loc[i, "seq"])
            this.__CreateDistanceFile(distance, splice_type)

    def distanceRI(this, splice_type: str = "RI_+") -> int:
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "RiStart", "RiEnd"]]  # récupère les coordonnées des sites de splicing
        distance = {
            "prot_start_splice_start": [],
            "prot_start_splice_end": [],
            "prot_end_splice_start": [],
            "prot_end_splice_end": [],
            "transcript_ID" : [], "prot_seq" : []
        }
        # Browse all the protein fixation sites
        for i in range(this.__data_prot.shape[0]):
            # Select the splicing sites of the same gene to avoid heavy computational loss
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            for j in range(splicing_same_gene.shape[0]):
                distance["prot_start_splice_start"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["RiStart"])
                distance["prot_start_splice_end"].append(this.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["RiEnd"])
                distance["prot_end_splice_start"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["RiStart"])
                distance["prot_end_splice_end"].append(this.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["RiEnd"])
                distance["transcript_ID"].append(this.__data_prot.loc[i, "ensembl_id"])
                distance["prot_seq"].append(this.__data_prot.loc[i, "seq"])
        this.__CreateDistanceFile(distance, splice_type)
        


    def ComputeDistance(this, splice_type : str = ""):
        
        
        distance_calculator = {
            "A5SS_+": this.distanceA5SS,
            "A5SS_-": this.distanceA3SS,
            "A3SS_+": this.distanceA3SS,
            "A3SS_-": this.distanceA5SS,
            "RI_-": this.distanceRI ,
            "RI_+": this.distanceRI #,
            # "MXE": this.distanceMXE,
            # "SE": this.distanceSE,
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
    splice_types = ["A5SS_+", "A5SS_-", "A3SS_+", "A3SS_-", "RI_+", "RI_-"]
    for splice_type in splice_types:
        dist.start(path_prot = "data_filteredFMRP.tsv", path_splicing = "filteredRmats", splice_type = splice_type)
        print(f"Distances for {splice_type} computed")