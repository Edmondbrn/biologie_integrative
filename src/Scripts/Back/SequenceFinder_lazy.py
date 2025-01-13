import pyensembl as pb
from pandas import read_csv, DataFrame
import os
import regex

NB_PROCESS = 4
ENSEMBL_NAME = "ensembl_id"
SEQUENCE_NAME = "seq"
RELEASE = 102

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class SequenceFinderLazy():
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

    def getDataProt(self):
        return self.__data_prot

    
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
            


        
    # On supprime la fonction d'alignement et on utilise directement start_ensembl / end_ensembl :
    def start(self):
        """
        Lance la conversion des positions splicées (start_ensembl, end_ensembl) en coordonnées génomiques
        sans réaliser d'alignement.
        """
        # Initialisation des listes où stocker les coordonnées génomiques calculées
        start_list = []
        end_list = []

        # Pour chaque ligne, on va chercher l'ID Ensembl, le start et end ensembl
        for i in range(len(self.__data_prot)):
            ensembl_id_val = self.__data_prot.loc[i, "ensembl_id"]
            start_rna = self.__data_prot.loc[i, "start_ensembl"]
            end_rna = self.__data_prot.loc[i, "end_ensembl"]

            try:
                # Récupération du transcript par ID
                transcript = self.__bdd.transcript_by_id(ensembl_id_val)
            except:
                # Si la récupération échoue, on enregistre None
                start_list.append("Not found")
                end_list.append("Not found")
                continue

            # Vérifie que les positions sont valides
            if isinstance(start_rna, int) and isinstance(end_rna, int):
                # On crée la plage de positions splicées
                spliced_positions = range(start_rna, end_rna + 1)
                # Conversion vers les positions génomiques
                genomic_range = self.__spliced_to_genomic(transcript, spliced_positions)

                if len(genomic_range) > 0:
                    # On prend le premier élément pour le début, et le dernier pour la fin
                    start_g = genomic_range[0][0]
                    end_g = genomic_range[-1][1]
                    start_list.append(start_g)
                    end_list.append(end_g)
                else:
                    start_list.append("Not found")
                    end_list.append("Not found")
            else:
                start_list.append("Not found")
                end_list.append("Not found")

        # On intègre les résultats dans le DataFrame
        self.__data_prot["start_genomic"] = start_list
        self.__data_prot["end_genomic"] = end_list

        # Sauvegarde du fichier
        # self.__data_prot.to_csv("data_filteredfinal3.tsv", sep="\t", index=False)
     
        




if __name__ == "__main__":
    df_prot = read_csv("data_filtered.tsv", sep = "\t", header = 0)
    app = SequenceFinder(df_prot)
    app.start()