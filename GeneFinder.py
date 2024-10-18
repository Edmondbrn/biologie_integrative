import re
import requests
from pyensembl import EnsemblRelease, Gene
from Bio import SeqIO
from io import StringIO
from os import chdir, path

# Changer le répertoire de travail
chdir(path.dirname(__file__))

class GeneFinder():
    """
    Class to convert RNA coordinates into DNA coordinates
    """

    def __init__(self, gene_name : str, RNA_seq : str, species : str ='human', release : int =104):
        self.__gene_name = gene_name
        self.__pattern = RNA_seq
        self.__species = species
        self.__release = release
        self.__data = EnsemblRelease(release, species = self.__species)
        self.__gene = None


    def __DownloadIndex(self):
        """Méthod to donwload and index the data if necessary"""
        self.__data.download()
        self.__data.index()
        self.__gene = self.__data.genes_by_name(self.__gene_name)

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

    def __DisplayGeneInfo(self):
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
        self.__response = requests.get(self.__url) # send a request to the API
        if self.__IsStatusOkay():
            self.__fasta_io = StringIO(self.__response.text) # create a StringIO object
            self.__record = SeqIO.read(self.__fasta_io, "fasta") # read the sequence in FASTA format
        
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
            print("Motif not found")
        else:
            if self.__IsUnique(self.__match): # if the motif is found only once
                self.__coordinates_start = self.__gene.start + self.__match[0].start()
                self.__coordinates_end = self.__gene.start + self.__match[0].end()
                print(f"Motif found at position {self.__coordinates_start} : {self.__coordinates_end}")
            else:
                print("Multiple matches found")

    def _Finder(self):
        """
        Method to find the gene, get the gene sequence and find a motif in the gene sequence
        """
        self.__DownloadIndex()
        self.__CheckGene()
        self.__GetGeneSequence()
        self.__DisplayGeneInfo()
        self.__FindMotif()


# Rechercher le gène Msh2
gene_name = "Dlc1"

mmotif = "CTGCAGATAGATTACAAGG"

finder = GeneFinder(gene_name, mmotif, "mouse")
finder._Finder()
