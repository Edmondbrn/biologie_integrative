from pyensembl import Genome

def get_cdna_and_exon_info(transcript_id, fasta_path, gtf_path, species='mouse'):
    # Initialiser l'objet Genome avec les chemins vers les fichiers FASTA et GTF
    genome = Genome(
        reference_name=species,
        annotation_name="custom_annotation",
        gtf_path_or_url=gtf_path,
        transcript_fasta_paths_or_urls=fasta_path,
        protein_fasta_paths_or_urls=None  # Si vous n'avez pas de fichier FASTA de protéines, mettez None
    )
    
    # Télécharger et indexer les données
    genome.download()
    genome.index()
    
    # Récupérer le transcrit à partir de son ID
    transcript = genome.transcript_by_id(transcript_id)
    print(transcript_id)
    # Récupérer la séquence cDNA
    cdna_sequence = genome.transcript_sequence(transcript_id)
    
    # Récupérer les informations sur les exons
    exons = transcript.exons
    
    return cdna_sequence, exons

# Exemple d'utilisation
fasta_path = "/home/edmond/aws/references/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa"  # Remplacez par le chemin vers votre fichier FASTA
gtf_path = "/home/edmond/aws/references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf"  # Remplacez par le chemin vers votre fichier GTF
transcript_id = "ENSMUST00000193812"  # Remplacez par votre identifiant d'ARNm Ensembl

cdna_sequence, exons = get_cdna_and_exon_info(transcript_id, fasta_path, gtf_path)

# Afficher la séquence cDNA
print(f"cDNA sequence for transcript {transcript_id}:\n{cdna_sequence[:100]}...")  # Afficher les 100 premiers nucléotides

# Afficher les informations des exons
print(f"\nExons for transcript {transcript_id}:")
for exon in exons:
    print(f"Exon {exon.exon_id}: start={exon.start}, end={exon.end}, length={exon.length}, sequence: {exon.sequence}")