import re
import requests
from pyensembl import EnsemblRelease
from Bio import SeqIO
from io import StringIO
from numba import jit
from os import chdir, path
from time import time

# Changer le répertoire de travail
chdir(path.dirname(__file__))

class GeneFinder():

    """
    Class to convert RNA coordinates into DNA coordinates
    """

    def __init__(self, gene_name : list[str], RNA_seq : list[str], species : str ='human', release : int = 113):
        if not isinstance(gene_name, list):
            raise ValueError("The gene name must be a list")
        if not isinstance(RNA_seq, list):
            raise ValueError("The RNA sequence must be a list")
        self.__gene_name = gene_name
        self.__pattern = RNA_seq
        self.__species = species
        self.__release = release
        self.__data = EnsemblRelease(self.__release, species = self.__species)
        self.__gene = None
        self.__session = requests.Session()

    def _destroy(self):
        """
        Methode to destroy the object
        """
        del self.__gene_name
        del self.__pattern
        del self.__species
        del self.__release
        del self.__data
        del self.__gene


    def _DownloadIndex(self):
        """Méthod to donwload and index the data if necessary"""
        self.__data.download()
        self.__data.index()

    def __IsNone(self, object):
        """
        Méthod to check if the gene is None
        """
        if not object:
            return True
        else:
            return False
        
    def __IsUnique(self, object):
        """
        Méthod to check if the gene is unique
        """
        if len(object) == 1:
            return True
        else:
            return False

    def __ChooseGene(self):
        """
        Methode to choose the gene if there are multiple matches
        """
        print(f"Multiple genes found for {self.__gene_name}.\n Please select one of the following gene IDs:")
        for index, gene in enumerate(self.__gene): # browse the list of genes
            print(f"N°{index}. Gene ID: {gene.gene_id}")
        try:
            self.__choice = int(input("Enter the number of the gene you want to select: "))
            if self.__choice not in range(len(self.__gene)):
                print("Invalid choice. Please enter a number corresponding to the gene you want to select.")
                self.__ChooseGene()
        except ValueError:
            print("Invalid choice. Please enter a number corresponding to the gene you want to select.")
            self.__ChooseGene()

    def __CheckGene(self):
        """
        Methode to check if the gene is in the database and if there are various matches
        """
        if self.__IsNone(self.__gene): # if there is no match
            print(f"Gene {self.__gene_name} not found in Ensembl release {self.__release}")
        else:
            if self.__IsUnique(self.__gene): # if there is only one match
                print(f"Gene {self.__gene_name} found in Ensembl release {self.__release}")
                self.__gene = self.__gene[0]
            else:
                self.__ChooseGene()
                self.__gene = self.__gene[self.__choice]

    def _DisplayGeneInfo(self):
        """
        Methode to display information about the gene
        """
        print(f"Name of the gene: {self.__gene.name}")
        print(f"Gene ID: {self.__gene.gene_id}")
        print(f"Chromosome: {self.__gene.contig}")
        print(f"Start position: {self.__gene.start}")
        print(f"End position: {self.__gene.end}")
        print(f"Strand: {self.__gene.strand}")

    def __IsStatusOkay(self):
        """
        Methode to check if the status code of the response is 200
        """
        if self.__response.status_code == 200:
            return True
        else:
            return False

    def __GetGeneSequence(self):
        """
        Method to get the gene sequence using the Ensembl REST API
        """
        self.__url = f"https://rest.ensembl.org/sequence/id/{self.__gene.gene_id}?content-type=text/x-fasta" # normalement renvoie le gène sur le brin codant donc l'ARN correspond à cet output mais à vérifier
        self.__response = self.__session.get(self.__url)  # Utiliser la session persistante
        if self.__IsStatusOkay():
            self.__fasta_io = StringIO(self.__response.text) # create a StringIO object
            self.__record = SeqIO.read(self.__fasta_io, "fasta") # read the sequence in FASTA format
    

    def __GetGeneSequencebis(self):
        self.__GENE = self.__data.gene_by_id(self.__gene.gene_id)
        print(self.__GENE)

    def __GetInfoGeneSequence(self):
        """
        Method to get information about the gene sequence
        """
        print(f"ID of the sequence: {self.__record.id}")
        print(f"Description: {self.__record.description}")
        print(f"Sequence: {self.__record.seq}")

    def __FindMotif(self):
        """
        Method to find a motif in the gene sequence
        """
        self.__match = list(re.finditer(self.__pattern, str(self.__record.seq)))
        if self.__IsNone(self.__match): # if the motif is not found
            return ("Pattern not found", "Pattern not found")
        else:
            if self.__IsUnique(self.__match): # if the motif is found only once
                self.__coordinates_start = self.__gene.start + self.__match[0].start()
                self.__coordinates_end = self.__gene.start + self.__match[0].end()
                return (self.__coordinates_start, self.__coordinates_end)
            else:
                return "Multiple matches found"
    
    def GetGeneEnsemblName(self):
        try:
            self.__gene = self.__data.genes_by_name(self.__gene_name)
            return True
        except:
            self.data_store[(self.__gene_name, self.__pattern, self.cpt)] = ("Gene not found", "Gene not found")
            self.cpt += 1
            return False
            
    def __BrowseGene(self):
        """
        Method to browse the list of genes and to search for the gene, get the gene sequence and find a motif in the gene sequence
        """
        self.data_store = dict()
        self.cpt = 0
        # self.length = len(self.__gene_name)
        for gene, motif in zip(self.__gene_name, self.__pattern):
            self.__gene_name = gene
            self.__pattern = motif
            self.start_1st_url = time()
            if not self.GetGeneEnsemblName() : # verify if ensembl succeed to find the gene
                continue
            self.end_1st_url = time()
            print(f"Gene {self.__gene_name} found in {self.end_1st_url - self.start_1st_url} seconds")
            self.__CheckGene() # check if the gene is in the database and if there are various matches
            self.start_2nd_url = time()
            self.__GetGeneSequence() # get the gene sequence using the Ensembl REST API
            self.end_2nd_url = time()
            print(f"Gene sequence of {self.__gene_name} downloaded in {self.end_2nd_url - self.start_2nd_url} seconds")
            self.start_sequencebis = time()
            self.__GetGeneSequencebis()
            self.end_sequencebis = time()
            print(f"LOCALLY Gene sequence of {self.__gene_name} downloaded in {self.end_sequencebis - self.start_sequencebis} seconds")
            self.start_searching = time()
            self.tuple_corrd = self.__FindMotif() # find the RNA sequence in the gene sequence
            self.end_searching = time()
            print(f"Motif {self.__pattern} found in {self.end_searching - self.start_searching} seconds")
            self.data_store[(self.__gene_name, self.__pattern, self.cpt)] = self.tuple_corrd
            self.cpt += 1
        #     print(f"Gene \r{self.cpt}/{self.length} processed", end = "")
        # print(f"Gene {self.length}/{self.length} processed")
        return self.data_store

   
    
    def __FormatResults(self):
        """
        Method to format the output file of the analysis
        """
        self.fh = open("results.txt", "w")
        self.fh.write("Gene\tMotif\tCoordinates_start\tCoordinates_end\n")
        for key, value in self.final_data.items():
            self.fh.write(f"{key[0]}\t{key[1]}\t{value[0]}\t{value[1]}\n")
        self.fh.close()
        print("Results saved in results.txt")




    def _Finder(self):
        """
        Method to find the gene, get the gene sequence and find a motif in the gene sequence
        """
        self._DownloadIndex()
        self.final_data = self.__BrowseGene()
        self.__FormatResults()

if __name__ == "__main__":
    # Rechercher le gène Msh2
    gene_name = ["Dlc1", "Msh2", "Hspa9"]
    mmotif = ["CTGCAGATAGATTACAAGG", "ATTTGGGTTAGCATGGGCT", "GCTCAAGGATTATGACTGC"]
    finder = GeneFinder(gene_name, mmotif, "mouse", 111)
    finder._Finder()
