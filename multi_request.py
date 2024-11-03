import requests
from  pyensembl import EnsemblRelease
from multiprocessing import Pool
import time
import json
import pandas as pd

GENE_PER_PROCESS = 3
NB_PROCESS = 16
REQUESTS_PER_SECOND = 15


 


def request_gene(parameters : tuple):
    trans_ids, start_list, end_list, specie = parameters
    # URL de base pour l'API REST d'Ensembl
    base_url = "https://rest.ensembl.org"
    dict_sequences = {}
    for transcrit_id, start, end in zip(trans_ids, start_list, end_list):
        # Construire l'URL pour la requête
        url = f"{base_url}/map/cdna/{transcrit_id}/{start}..{end}"
        # Envoyer la requête GET
        response = requests.get(url, headers={"Content-Type" : "application/json"})
        
        # Vérifier le statut de la réponse
        if response.status_code == 200:
            # TODO Voir pour diviser le cas où on a plusieurs fois le même transcrit
            dict_sequences[transcrit_id] = response.text.strip()
        else:
            print(f"Failed to retrieve transcript ID: {response.status_code}")
            dict_sequences[transcrit_id] = None
    return dict_sequences
            
def get_dna_sequences(transcript_ids_list : list[str] , start_trans :list[int], end_trans : list[int], specy : str = "mouse"):

    id = 0
    nb_id = len(transcript_ids_list)
    params = list()
    while id < nb_id : # On découpe la liste des ensembl_ids en batchs de taille NB_IDS_PER_REQUEST
        params.append((transcript_ids_list[id:id+GENE_PER_PROCESS], start_trans[id:id+GENE_PER_PROCESS], end_trans[id:id+GENE_PER_PROCESS], specy))
        id += GENE_PER_PROCESS

    pool = Pool(processes = NB_PROCESS) # On crée un pool de processus pour la parallélisation
    gene_sequences : list[dict] = pool.map(request_gene, params) # on exécute les requêtes en parallèle
    final_dict = dict()
    for dict_genes  in gene_sequences:
        final_dict = {**final_dict, **dict_genes}
    return final_dict


# Exemple d'utilisation
transcript_ids_list = [
    "ENSMUST00000193812", "ENSMUST00000193813", "ENSMUST00000193814", "ENSMUST00000193815", "ENSMUST00000193816",
    "ENSMUST00000193817", "ENSMUST00000193818", "ENSMUST00000193819", "ENSMUST00000193820", "ENSMUST00000193821",
    "ENSMUST00000193822", "ENSMUST00000193823", "ENSMUST00000193824", "ENSMUST00000193825", "ENSMUST00000193826",
    "ENSMUST00000193827", "ENSMUST00000193828", "ENSMUST00000193829", "ENSMUST00000193830", "ENSMUST00000193831"
]

start_trans = [
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    1, 1, 1, 1, 1
]

end_trans = [
    1000, 1100, 1200, 1300, 1400,
    1500, 1600, 1700, 1800, 1900,
    2000, 2100, 2200, 2300, 2400,
    2500, 2600, 2700, 2800, 2900
]
transcript_dict = {transcript_id: (start, end) for transcript_id, start, end in zip(transcript_ids_list, start_trans, end_trans)}

# Appeler la fonction get_dna_sequences
sequences = get_dna_sequences(transcript_ids_list, start_trans, end_trans, specy="mouse")

# Afficher les séquences récupérées
for transcript_id, sequence in sequences.items():
    if sequence != None :
        sequence = json.loads(sequence) # convertit la str en dictionnaire
        start_gen = sequence["mappings"][0]["start"]
        start_arn = transcript_dict[transcript_id][0]
        end_gen = sequence["mappings"][0]["end"]
        end_arn = transcript_dict[transcript_id][1]
        print(f"Transcript ID: {transcript_id} position {start_arn} : {end_arn}, Genomic position: {start_gen} : {end_gen}")  # Afficher les 100 premiers nucléotides
    else:
        print(f"Transcript ID: {transcript_id} not found")
