import os
import pandas as pd

os.chdir(os.path.dirname(__file__))

class Distances():

    """
    This class is for computing distances between splicing sites and protein fixation site on mRNA
    """


    def __init__(self):
        self.__data_prot : pd.DataFrame = None
        self.__data_splicing : dict[ str : pd.DataFrame ] = None
    

    def __LoadDataProt(self, file_name : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data of the protein
        file_name : name of the file of the protein
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        self.__data_prot = pd.read_csv(file_name, sep = sep)
        return self.__data_prot
    
    def __LoadDataSplicing(self, path : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data from alternative splicing
        path : path to the folder containing the files
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        for file in os.listdir(path):
            self.__data_splicing[file] = pd.read_csv(f"{path}/{file}", sep = sep)
    
    distance_calculator = {
        "A5SS_+": distanceA5SS,
        "A5SS_-": distanceA3SS,
        "A3SS_+": distanceA3SS,
        "A3SS_-": distancerA5SS,
        "RI": distanceRI,
        "MXE": distanceMXE,
        "SE": distanceSE,
        }

    def distanceA5SS(self):
        splicing_coordinates = self.__data_splicing["A5SS_+"][["shortSplice", "longSplice", "shareSplice"]] #récupère les coordonnées des sites de splicing
        protein_coordinates = self.__data_prot #récupère les coordonnées du site de fixation
        

    def ComputeDistance(self):
        





if __name__ == "__main__":
    dist = Distances()
    