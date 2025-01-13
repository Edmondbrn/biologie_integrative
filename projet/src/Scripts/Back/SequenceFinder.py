import pyensembl as pb
from pandas import read_csv, DataFrame
import os
import regex
import ast

NB_PROCESS = 4
ENSEMBL_NAME = "ensembl_id"
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
        litteral_col = []
        for row in self.__data_prot["ensembl_id"]:
            litteral_col.append(ast.literal_eval(row))
        self.__data_prot["ensembl_id"] = litteral_col
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
    
    @staticmethod
    def isRnaCoordNumber(coord):
        return isinstance(coord, int)
    

    


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

                dict_tmp[(tuple(ensembl_id[j]), j)] = fixation_sequence[j]
                # On a un dictionnaire avec en clé les identifiants Ensembl et en valeur les séquences de fixation
            self.__general_list.append(dict_tmp)
        return None
    
    @staticmethod
    def align(cDNA : str, sequence : str) -> tuple[int, int]:
        """
        Method to align the sequence on the cDNA
        """
        pattern = f"({sequence}){{e<=0}}"
        matches = list(regex.finditer(pattern, cDNA))
        
        if not matches:
            return ("Not found", "Not found")
        elif len(matches) == 1:
            start = matches[0].start()
        else:
            print("Multiple start found")
            for idx, match in enumerate(matches):
                print(f"{idx}: {match.start()}")
            position = int(input(f"Choose which one you want to keep: "))
            start = matches[position].start()
        
        end = start + len(sequence)
        return start, end

    def __alignSequences(self, parameters: list[pb.Database, dict[tuple[str, int], str]]) -> dict[tuple[str, int], tuple[int, int]]:
        """
        Method to align the sequences on the cDNA from pyensembl
        """
        bdd, input_data = parameters
        ensembl_ids = input_data.keys()
        sequences = input_data.values()
        result = {}
        
        for ensembl_id, sequence in zip(ensembl_ids, sequences):
            transcript = None
            for ensembl_id_part in ensembl_id[0]:
                try:
                    transcript = bdd.transcript_by_id(ensembl_id_part)
                    break  # Si on trouve un transcript valide, on sort de la boucle
                except:
                    continue  # Si on ne trouve pas, on continue avec le prochain identifiant
            
            if transcript is None:
                result[ensembl_id] = ("unknown", "unknown"), [("unknown", "unknown")], 'unknown'
            else:
                cDNA = transcript.sequence
                gene = transcript.gene_id
                rna_start, rna_end = SequenceFinder.align(cDNA, sequence)
                if SequenceFinder.isRnaCoordNumber(rna_start) and SequenceFinder.isRnaCoordNumber(rna_end):
                    self.genomic_coordinate_list = self.__spliced_to_genomic(transcript, range(rna_start, rna_end + 1))
                    result[ensembl_id] = ((rna_start, rna_end), self.genomic_coordinate_list, gene)
                else:
                    result[ensembl_id] = ("Not found", "Not found"), [("Not found", "Not found")], "Not found"
        
        return result
    
    def __spliced_to_genomic(self, transcript : pb.Transcript, spliced_positions : list[int]):
        """
        Convert a spliced position to a genomic position.
        take the concerned transcript as argument
        and a list of spliced position
        """
        current_spliced_position = 0
        exon_number = 0
        exons = list(transcript.exons)
        genomic_positions = []
     
        for spliced_position in spliced_positions: # browse the list of spliced position
            if spliced_position == "Not found":
                genomic_positions.append(None)
                continue
            for i in range(exon_number, len(exons)): # browse the exon but we kept the exon number to avoid to browse all the exons each time
                exon = exons[i]
                exon_length = exon.end - exon.start
                if current_spliced_position + exon_length > spliced_position:
                    offset = spliced_position - current_spliced_position
                    genomic_positions.append(exon.start + offset)
                    break
                else:
                    current_spliced_position += exon_length
                    exon_number += 1
                if exon_number > len(exons):
                    genomic_positions.append(None)  # if the spliced position is out of range
        genomic_range = self.__determineGap(genomic_positions)
        return genomic_range

    def __determineGap(self, genomic_positions : list[int]) -> list[int]:
        """
        Method to determine the gap between genomic positions thus to determine 
        if the fixation sequences is between various exons or not
        """
        min_pos = 0
        self.__genomic_range_list : list[tuple[int, int]] = []
        for j in range(len(genomic_positions) - 1): # browse all the genomic positions
            if genomic_positions[j] is None: # if we have None position we kept the information
                self.__genomic_range_list.append((None, None))
            elif genomic_positions[j + 1] is None:
                self.__genomic_range_list.append((genomic_positions[j], None))
            else:
                gap = genomic_positions[j + 1] - genomic_positions[j]
            if gap == 1 and min_pos == 0:
                min_pos = genomic_positions[j]
            elif gap > 1 and min_pos != 0 or j == len(genomic_positions) - 2: # si on a un gap de plus de 1 ou si on est à la fin de la liste contenant les positions
                max_pos = genomic_positions[j+1] if j == len(genomic_positions) - 2 else genomic_positions[j]
                self.__genomic_range_list.append((min_pos, max_pos))
                min_pos = 0
        return self.__genomic_range_list
            
    def __addRnaCoordinates(self, coord_dict : dict) -> None:
        """
        Method to add the coordinates to the DataFrame
        """
        start_list = {"start_ensembl" : list()}
        end_list = {"end_ensembl" : list()}
        gene_list_id = {"GeneID" : list()}
        for id, coord in coord_dict.items():
            start_list["start_ensembl"].append(coord[0][0])
            end_list["end_ensembl"].append(coord[0][1])
            gene_list_id["GeneID"].append(coord[2])

        self.__data_prot = self.__data_prot.join(DataFrame(start_list))
        self.__data_prot = self.__data_prot.join(DataFrame(end_list))
        self.__data_prot = self.__data_prot.join(DataFrame(gene_list_id))
        return None
    
    def __addGenomicCoordinates(self, coord_list : dict) -> None:
        """
        Method to add the genomic coordinates to the DataFrame
        """
        start_list = {"start_genomic_complete" : list()}
        end_list = {"end_genomic_complete" : list()}
        for coord in coord_list.keys():
            start_tuple = list()
            end_tuple = list()
            for tuple_position in coord_list[coord][1]:
                start_tuple.append(tuple_position[0])
                end_tuple.append(tuple_position[1])
            start_list["start_genomic_complete"].append(start_tuple)
            end_list["end_genomic_complete"].append(end_tuple)
        self.__data_prot = self.__data_prot.join(DataFrame(start_list))
        self.__data_prot = self.__data_prot.join(DataFrame(end_list))
        return None

    
    def start(self):
        self.___getFixationSequence()
        # pool = Pool(NB_PROCESS)
        # réunis les dictionnaires séparés dans la liste générale (artefact du multithreading si jamais on veut en faire)
        all_dict = {}
        for element in self.__general_list:
            all_dict.update(element)

        __dict_coord = self.__alignSequences((self.__bdd, all_dict))
        self.__addRnaCoordinates(__dict_coord)
        self.__addGenomicCoordinates(__dict_coord)
        start_list = {"start_genomic" : list()}
        end_list = {"end_genomic" : list()}
        for i in range(len(self.__data_prot)):
            start_list["start_genomic"].append(self.__data_prot.iloc[i]["start_genomic_complete"][0])
            end_list["end_genomic"].append(self.__data_prot.iloc[i]["end_genomic_complete"][-1])
        self.__data_prot = self.__data_prot.join(DataFrame(start_list))
        self.__data_prot = self.__data_prot.join(DataFrame(end_list))
        self.__data_prot.to_csv("data_filteredfinal2.tsv", sep = "\t", index = False)

        
     
        




if __name__ == "__main__":
    df_prot = read_csv("/home/edmond/Documents/GB5/biologie_integrative/projet/src/Ressources/data/FMRP_Binding_sites_mouse_Maurin_NAR_2014_merged.tsv_converted", sep = "\t", header = 0)
    app = SequenceFinder(df_prot)
    app.start()