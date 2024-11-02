from Bio import SeqIO
import time
import psutil

def print_memory_usage():
    process = psutil.Process()
    mem_info = process.memory_info()
    print(f"Memory usage: {mem_info.rss / 1024 ** 2:.2f} MB")

# Afficher l'utilisation de la mémoire avant l'indexation
print("Memory usage before indexing:")
print_memory_usage()

# Indexer le fichier FASTA
record_dict = SeqIO.index("/home/edmond/aws/references/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa", "fasta")

# Afficher l'utilisation de la mémoire après l'indexation
print("Memory usage after indexing:")
print_memory_usage()

def get_exon_sequence(chromosome, start, end):
    # Récupérer la séquence du chromosome
    chromosome_seq = record_dict[chromosome]
    
    # Extraire la séquence de l'exon
    exon_sequence = chromosome_seq.seq[start-1:end]  # Les indices sont 0-based, donc on soustrait 1 de start
    return exon_sequence

# Exemple d'utilisation
chromosome = "1"  # Remplacez par le chromosome de votre exon
start = 1000000  # Remplacez par la position de début de votre exon
end = 1001000  # Remplacez par la position de fin de votre exon

debut = time.time()
exon_sequence = get_exon_sequence(chromosome, start, end)
print(f"Exon sequence: {exon_sequence}")
print(time.time() - debut)