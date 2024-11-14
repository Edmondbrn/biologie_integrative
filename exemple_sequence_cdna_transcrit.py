

from pyensembl import EnsemblRelease

def spliced_to_genomic(transcript, spliced_position):
    """
    Convert a spliced position to a genomic position.
    """
    current_spliced_position = 0
    for exon in list(transcript.exons):
        exon_length = exon.end - exon.start + 1
        if current_spliced_position + exon_length > spliced_position:
            offset = spliced_position - current_spliced_position
            return exon.start + offset
        current_spliced_position += exon_length
    raise ValueError("Spliced position is out of range")

def get_cdna_and_exon_info(transcript_id, species='mouse'):
    # Initialiser pyensembl pour la souris (Mus musculus)
    data = EnsemblRelease(102, species=species)
    data.download()
    data.index()
    
    # Récupérer le transcrit à partir de son ID
    transcript = data.transcript_by_id(transcript_id)
    print(transcript.exon_intervals)
    
    # Exemple de conversion d'une position chromosomique en position dans l'ARNm épissé
    chromosomal_position = 13574523
    spliced_position = transcript.spliced_offset(chromosomal_position)
    print(f"Spliced position for chromosomal position {chromosomal_position}: {spliced_position}")
    
    # Exemple de conversion d'une position dans l'ARNm épissé en position chromosomique
    spliced_position = 477
    genomic_position = spliced_to_genomic(transcript, spliced_position)
    print(f"Genomic position for spliced position {spliced_position}: {genomic_position}")
    
    # Afficher les 100 premiers nucléotides de la séquence cDNA
    print(transcript.sequence[:100])

# Exemple d'utilisation
trans_id = "ENSMUST00000191615"  # Remplacez par votre identifiant de transcrit Ensembl
get_cdna_and_exon_info(trans_id)