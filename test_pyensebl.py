import requests
from  pyensembl import EnsemblRelease
from multiprocessing import Pool
import time


def request_gene(parameters : tuple):
    gene_ids, specie, data_release, nb_flancant = parameters
    # URL de base pour l'API REST d'Ensembl
    base_url = "https://rest.ensembl.org"
    dict_sequences = {}
    for gene_id in gene_ids:
        gene = data_release.gene_by_id(gene_id)
        # Construire l'URL pour la requête
        url = f"{base_url}/sequence/region/{specie}/{gene.contig}:{gene.start - nb_flancant}..{gene.end + nb_flancant}:{gene.strand}?content-type=text/plain"
        # Envoyer la requête GET
        response = requests.get(url)
        
        # Vérifier le statut de la réponse
        if response.status_code == 200:
            dict_sequences[gene_id] = response.text.strip()
        else:
            print(f"Failed to retrieve sequence: {response.status_code}, {response.text}")
            dict_sequences[gene_id] = None
        return dict_sequences
            

def get_dna_sequences(gene_ids_list : list[str], nb_flancant : int = 600 , specy : str = "human"):

    id = 0
    nb_id = len(gene_ids_list)
    params = list()
    data_release = EnsemblRelease(species = specy)

    while id < nb_id : # On découpe la liste des ensembl_ids en batchs de taille NB_IDS_PER_REQUEST
        params.append((gene_ids_list[id:id+GENE_PER_PROCESS], specy, data_release, nb_flancant))
        id += GENE_PER_PROCESS

    pool = Pool(processes = NB_PROCESS) # On crée un pool de processus pour la parallélisation
    gene_sequences : list[dict] = pool.map(request_gene, params) # on exécute les requêtes en parallèle
    
    final_dict = dict()
    for dict_genes  in gene_sequences:
        final_dict = {**final_dict, **dict_genes}
    return final_dict


# Exemple d'utilisation
ensembl_ids = [
    "ENSMUSG00000064372", "ENSMUSG00000064371", "ENSMUSG00000064370", "ENSMUSG00000064369", "ENSMUSG00000064368",
    "ENSMUSG00000064367", "ENSMUSG00000064366", "ENSMUSG00000064365", "ENSMUSG00000064364", "ENSMUSG00000064363"
]
data = EnsemblRelease(species = "mouse")
dict_gene_info = {}
for ensembl_d in ensembl_ids:
    dict_gene_info[ensembl_d] = data.gene_by_id(ensembl_d)
NT_FLANCANT = 600
GENE_PER_PROCESS = 1
NB_PROCESS = 8

# gene par gene 

# debut = time.time()
# for gene in dict_gene_info.values():
#     sequence = get_dna_sequences(gene.contig, gene.start -600, gene.end +600, gene.strand, "mouse")
#     # if sequence:
#     #     print(f"ID : {gene.gene_id} Sequence: {sequence[:100]}...")  # Afficher les 100 premiers nucléotides
# print(time.time() - debut)

# Multiprocess
debut = time.time()
print(get_dna_sequences(ensembl_ids, nb_flancant= 600, specy="mouse"))
print(time.time() - debut)