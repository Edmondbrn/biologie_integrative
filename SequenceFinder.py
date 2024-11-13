import pyensembl as pb
from pandas import read_csv, DataFrame
import os
import re
from multiprocessing import Pool
import time

NB_PROCESS = 4
ENSEMBL_NAME = "ensembl_transcript_id"
SEQUENCE_NAME = "seq"
RELEASE = 102

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class SequenceFinder():
    """
    Cette classe va permettre de récupérer les coordonnées ARN données en entrée et de les convertir en 
    coordonnées ADN à l'aide de l'API REST de Ensembl.
    """


    def __init__(self, 
                #  data_splice : DataFrame, 
                 data_prot : DataFrame, 
                 aim_assembly : str = "GRCm38", 
                 species : str = "mouse"):
        """
        Constructeur de la classe SequenceFinder.
        :param data_splice: DataFrame contenant les coordonnées ARN à convertir.
        :param data_prot: DataFrame contenant les coordonnées protéiques à convertir.
        :param species: Espèce sur laquelle effectuer la recherche.
        """
        # self.__data_splice = data_splice
        self.__data_prot = data_prot
        self.__species = species
        self.__assembly = aim_assembly
        # TODO change the ensemble release according to user input
        self.__bdd = pb.EnsemblRelease(RELEASE, species=self.__species) # version 102 pour avoir l'assemblage 38 de la souris
        self.__bdd.download()
        self.__bdd.index()

    
    def setAttribute(self, attibute : str, value):
        """
        Cette méthode va permettre de définir les attributs de la classe.
        """
        setattr(self, attibute, value)
        return None

    def getAttribute(self, attribute : str):
        """
        Cette méthode va permettre de récupérer les attributs de la classe.
        """
        return getattr(self, attribute)
    
    @staticmethod
    def isUnique(element):
        return len(element) == 1
    
    @staticmethod
    def isNone(element):
        return len(element) == 0
    

    


    def ___getFixationSequence(self) -> list[ dict[ tuple[str, int] : str] ]:
        """
        Method to get the sequences and Ensembl ID from the file stored in self.data_prot
        Return a list of dictionnaries ti later use it in multiprocessing algorithm:
        Key : tuple with ensembl_ID and an index (we can have sevral times the same ensembl ID)
        Value : sequence
        """
        ensembl_id, fixation_sequence = self.__data_prot[ENSEMBL_NAME].to_list(), self.__data_prot[SEQUENCE_NAME].to_list()
        self.__general_list = []
        size_ensembl_id = len(ensembl_id)
        number_id_per_dict = size_ensembl_id // NB_PROCESS + NB_PROCESS
        for i in range(NB_PROCESS): # creation of the dictionnaries for the multiprocessing
            dict_tmp = {}
            for j in range(i * number_id_per_dict, (i + 1) * number_id_per_dict): # filling the dictionnaries
                if j >= size_ensembl_id: # if we reach the end of the list of ensembl ID (avoid out of range if the division is not absolute)
                    break
                dict_tmp[(ensembl_id[j], j)] = fixation_sequence[j]
            self.__general_list.append(dict_tmp)
        return None
    
    @staticmethod
    def align(cDNA : str, sequence : str) -> tuple[int, int]:
        """
        Method to align the sequence on the cDNA
        """
        matches = list(re.finditer(sequence, cDNA))
        
        if SequenceFinder.isNone(matches):
            return ("Not found", "Not found")
        elif SequenceFinder.isUnique(matches):
            start = matches[0].start()
        else:
            print("Multiple start found")
            for idx, match in enumerate(matches):
                print(f"{idx}: {match.start()}")
            position = int(input(f"Choose which one you want to keep: "))
            start = matches[position].start()
        
        end = start + len(sequence)
        return start, end
    
    @staticmethod
    def alignSequences(parameters : list[pb.Database , dict[ tuple[str, int] : str]]) -> dict[ tuple[str, int] : tuple[int, int]]:
        """
        Method to align the sequences on the cDNA from pyensembl
        """
        bdd, input_data = parameters
        ensembl_id = input_data.keys()
        sequences = input_data.values()
        result = {}
        for ensembl_id, sequence in zip(ensembl_id, sequences):
            try:
                transcript = bdd.transcript_by_id(ensembl_id[0])
            except:
                result[ensembl_id] = ("unknown_id", "unknown_id")
            else:
                cDNA = transcript.sequence
                start, end = SequenceFinder.align(cDNA, sequence)
                result[ensembl_id] = (start, end)
        return result
    
    
    def start(self):
        self.___getFixationSequence()
        # pool = Pool(NB_PROCESS)
        # réunis les dictionnaires séparés dans la liste générale

        all_dict = {}
        for element in self.__general_list:
            all_dict.update(element)

        liste = SequenceFinder.alignSequences((self.__bdd, all_dict))
        start_list = {"start_ensembl" : list()}
        end_list = {"end_ensembl" : list()}
        for coord in liste.values():
            start_list["start_ensembl"].append(coord[0])
            end_list["end_ensembl"].append(coord[1])

        self.__data_prot = self.__data_prot.join(DataFrame(start_list))
        self.__data_prot = self.__data_prot.join(DataFrame(end_list))
        self.__data_prot.to_csv("data_filtered2.tsv", sep = "\t", index = False)
     
        





df_prot = read_csv("data_filtered.tsv", sep = "\t", header = 0)
app = SequenceFinder(df_prot)
app.start()