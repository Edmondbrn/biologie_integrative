import os
import pandas as pd
import ast

os.chdir(os.path.dirname(__file__))

class Distances():

    """
    This class is for computing distances between splicing sites and protein fixation site on mRNA
    """


    def __init__(self):
        self.__data_prot : pd.DataFrame = pd.DataFrame()
        self.__data_splicing : dict[ str : pd.DataFrame ] = dict()
    

    def __LoadDataProt(self, file_name : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data of the protein
        file_name : name of the file of the protein
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        self.__data_prot = pd.read_csv(file_name, sep = sep)
        self.__data_prot["start_genomic"] = self.__data_prot["start_genomic"].apply(ast.literal_eval)
        self.__data_prot["end_genomic"] = self.__data_prot["end_genomic"].apply(ast.literal_eval)
        self.__FilterDataProt()

    def __FilterDataProt(self):
        # enlever les lignes avec des str sur la colonne start_ensembl
        self.__data_prot = self.__data_prot[(self.__data_prot["start_ensembl"] != "Not found") | (self.__data_prot["start_ensembl"] != "unknown")]

    
    def __LoadDataSplicing(self, path : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data from alternative splicing
        path : path to the folder containing the files
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        for file in os.listdir(path):
            self.__data_splicing[file.split(".")[0]] = pd.read_csv(f"{path}/{file}", sep = sep)
    
    
    def distanceA5SS(self, splice_type : str = "A5SS_+"):
        splicing_coordinates : pd.DataFrame = self.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]] #récupère les coordonnées des sites de splicing
        # Browse all the protein fixation sites
        for i in range(self.__data_prot.shape[0]):
            # Select the splicing sites of the same gene to avoid heavy computationnal loss
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == self.__data_prot.loc[i, "GeneID"]]
            # Browse all the splicing sites of the same gene
            distance = {"prot_start_short_splice_start" : [], "prot_start_short_splice_end" : [], "prot_end_short_splice_start" : [], "prot_end_short_splice_end" : [],
                        "prot_start_long_splice_start" : [], "prot_start_long_splice_end" : [], "prot_end_long_splice_start" : [], "prot_end_long_splice_end" : []}
            for j in range(splicing_same_gene.shape[0]):
                # Compute the distance between the protein fixation site and the splicing site
                distance["prot_start_short_splice_start"].append(self.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.loc[j, "shortSplice"])
                distance["prot_start_short_splice_end"].append(self.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.loc[j, "sharesplice"])
                distance["prot_end_short_splice_start"].append(self.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.loc[j, "shortSplice"])
                distance["prot_end_short_splice_end"].append(self.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.loc[j, "sharesplice"])
                distance["prot_start_long_splice_start"].append(self.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.loc[j, "longSplice"])
                distance["prot_start_long_splice_end"].append(self.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.loc[j, "sharesplice"])
                distance["prot_end_long_splice_start"].append(self.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.loc[j, "longSplice"])
                distance["prot_end_long_splice_end"].append(self.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.loc[j, "sharesplice"])

                
                
    def distanceA3SS(self, splice_type : str = "A3SS_+") :
            splicing_coordinates : pd.DataFrame = self.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]] #récupère les coordonnées des sites de splicing
            # Browse all the protein fixation sites
            for i in range(self.__data_prot.shape[0]):
                # Select the splicing sites of the same gene to avoid heavy computationnal loss
                splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == self.__data_prot.loc[i, "GeneID"]]
                # Browse all the splicing sites of the same gene
                distance = {"prot_start_short_splice_start" : [], "prot_start_short_splice_end" : [], "prot_end_short_splice_start" : [], "prot_end_short_splice_end" : [],
                        "prot_start_long_splice_start" : [], "prot_start_long_splice_end" : [], "prot_end_long_splice_start" : [], "prot_end_long_splice_end" : []}
                
                for j in range(splicing_same_gene.shape[0]):
                # Compute the distance between the protein fixation site and the splicing site
                # "shortSplice", "longSplice", "shareSplice"
                # Compute the distance between the protein fixation site and the splicing site
                    distance["prot_start_short_splice_start"].append(self.__data_prot.loc[i, "start_genomic"] - splicing_same_gene.loc[j, "shareSplice"])
                    distance["prot_start_short_splice_end"].append(self.__data_prot.loc[i, "start_genomic"] - splicing_same_gene.loc[j, "shortSplice"])
                    distance["prot_end_short_splice_start"].append(self.__data_prot.loc[i, "end_genomic"] - splicing_same_gene.loc[j, "shareSplice"])
                    distance["prot_end_short_splice_end"].append(self.__data_prot.loc[i, "end_genomic"] - splicing_same_gene.loc[j, "shortSplice"])
                    distance["prot_start_long_splice_start"].append(self.__data_prot.loc[i, "start_genomic"] - splicing_same_gene.loc[j, "shareSplice"])
                    distance["prot_start_long_splice_end"].append(self.__data_prot.loc[i, "start_genomic"] - splicing_same_gene.loc[j, "longSplice"])
                    distance["prot_end_long_splice_start"].append(self.__data_prot.loc[i, "end_genomic"] - splicing_same_gene.loc[j, "shareSplice"])
                    distance["prot_end_long_splice_end"].append(self.__data_prot.loc[i, "end_genomic"] - splicing_same_gene.loc[j, "longSplice"])


    def distanceRI(self, splice_type: str = "RI_+") -> int:
        splicing_coordinates: pd.DataFrame = self.__data_splicing[splice_type][["GeneID", "chr", "strand", "RiStart", "RiEnd"]]  # récupère les coordonnées des sites de splicing
        distance = {
            "prot_start_splice_start": [],
            "prot_start_splice_end": [],
            "prot_end_splice_start": [],
            "prot_end_splice_end": []
        }
            # Browse all the protein fixation sites
        print(self.__data_prot.shape[0])
        for i in range(self.__data_prot.shape[0]):
            # Select the splicing sites of the same gene to avoid heavy computational loss
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == self.__data_prot.loc[i, "GeneID"]]
            
            for j in range(splicing_same_gene.shape[0]):
                print(self.__data_prot.loc[i, "start_genomic"][0])
                print(splicing_same_gene.iloc[j]["RiStart"])
                distance["prot_start_splice_start"].append(self.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["RiStart"])
                distance["prot_start_splice_end"].append(self.__data_prot.loc[i, "start_genomic"][0] - splicing_same_gene.iloc[j]["RiEnd"])
                distance["prot_end_splice_start"].append(self.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["RiStart"])
                distance["prot_end_splice_end"].append(self.__data_prot.loc[i, "end_genomic"][-1] - splicing_same_gene.iloc[j]["RiEnd"])
        
        with open(f"distance{splice_type}.txt", "w") as f:
            for key in distance.keys():
                f.write(f"{key} : {distance[key]}\n")

        return distance
        
    

    def ComputeDistance(self, splice_type : str = ""):
        
        
        distance_calculator = {
            "A5SS_+": self.distanceA5SS,
            "A5SS_-": self.distanceA3SS,
            "A3SS_+": self.distanceA3SS,
            "A3SS_-": self.distanceA5SS,
            "RI_-": self.distanceRI ,
            "RI_+": self.distanceRI #,
            # "MXE": self.distanceMXE,
            # "SE": self.distanceSE,
            }

        return distance_calculator[splice_type](splice_type)
        
    def start(self, path_prot : str, path_splicing : str, splice_type : str):
        """
        Method to start the computation of the distances
        path_prot : path to the file containing the protein data with the csv file at the end
        path_splicing : path to the folder containing the splicing data
        """
        self.__LoadDataProt(path_prot)
        self.__LoadDataSplicing(path_splicing)
        print(self.ComputeDistance(splice_type))

if __name__ == "__main__":
    dist = Distances()
    dist.start(path_prot = "data_filteredFMRP.tsv", path_splicing = "filteredRmats", splice_type = "RI_-")
