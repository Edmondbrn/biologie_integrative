from pyensembl import EnsemblRelease

def get_cdna_and_exon_info(gene_id, species='mouse'):
    # Initialiser pyensembl pour la souris (Mus musculus)
    data = EnsemblRelease(104, species=species)
    data.download()
    data.index()
    
    # Récupérer le gène à partir de son ID
    gene = data.gene_by_id(gene_id)
    
    # Récupérer les transcrits associés au gène
    transcripts = gene.transcripts
    
    # Récupérer la séquence cDNA et les informations sur les exons pour chaque transcrit
    cdna_info = {}
    for transcript in transcripts:
        cdna_sequence = transcript.sequence
        exons = transcript.exons
        cdna_info[transcript.transcript_id] = {
            'cdna_sequence': cdna_sequence,
            'exons': exons
        }
    
    return cdna_info

# Exemple d'utilisation
gene_id = "ENSMUSG00000068457"  # Remplacez par votre identifiant de gène Ensembl

cdna_info = get_cdna_and_exon_info(gene_id)

# Afficher les informations pour chaque transcrit
for transcript_id, info in cdna_info.items():
    cdna_sequence = info['cdna_sequence']
    exons = info['exons']
    
    # Afficher la séquence cDNA
    print(f"cDNA sequence for transcript {transcript_id}:\n{cdna_sequence[:100]}...")  # Afficher les 100 premiers nucléotides
    
    # Afficher les informations des exons
    print(f"\nExons for transcript {transcript_id}:")
    for exon in exons:
        print(f"Exon {exon.exon_id}: start={exon.start}, end={exon.end}, length={exon.length}")