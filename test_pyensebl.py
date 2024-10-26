import requests
from  pyensembl import EnsemblRelease
from multiprocessing import Pool
import time
import pandas as pd

NT_FLANCANT = 600
GENE_PER_PROCESS = 10
NB_PROCESS = 8
REQUESTS_PER_SECOND = 15


def request_gene(parameters : tuple):
    gene_ids, specie, data_release, nb_flancant = parameters
    # URL de base pour l'API REST d'Ensembl
    base_url = "https://rest.ensembl.org"
    dict_sequences = {}
    for gene_id in gene_ids:
        if gene_id in dict_sequences.keys(): # saute si on a déjà récupéré la séquence d'un gène précis
            continue
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
        # Pause pour respecter la limite de 15 requêtes par seconde
        time.sleep(0.3)
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
df = pd.read_csv("/home/edmond/Documents/GB5/biologie_integrative/rmats_post/A5SS.MATS.JCEC.txt", 
                 sep = "\t")
ensembl_ids = df["GeneID"][:50].values
print(f"Il y a {len(pd.unique(ensembl_ids))} séquences uniques")
data = EnsemblRelease(species = "mouse")
dict_gene_info = {}
gene_to_remove = []
for ensembl_d in ensembl_ids:
    try :
        dict_gene_info[ensembl_d] = data.gene_by_id(ensembl_d)
    except:
        gene_to_remove.append(ensembl_d)
        continue
for gene in gene_to_remove:
    ensembl_ids.remove(gene)

# gene par gene 

# debut = time.time()
# for gene in dict_gene_info.values():
#     sequence = get_dna_sequences(gene.contig, gene.start -600, gene.end +600, gene.strand, "mouse")
#     # if sequence:
#     #     print(f"ID : {gene.gene_id} Sequence: {sequence[:100]}...")  # Afficher les 100 premiers nucléotides
# print(time.time() - debut)

# Multiprocess
debut = time.time()
dico = get_dna_sequences(ensembl_ids, nb_flancant= 600, specy="mouse")
print(f"Pour {len(dico)} sequences, le programme a mis {time.time() - debut} secondes")

fh = open("gene_sequences_new.txt", "w")
for id, sequences in dico.items():
    fh.write(f"{id}\n")
    try :
        for nuc in range(0, len(sequences) -len(sequences)%80, 80):
            fh.write(f'{sequences[0+nuc : nuc +80].strip()}\n')
    except:
        fh.write(f"{sequences}\n")


fh.close()