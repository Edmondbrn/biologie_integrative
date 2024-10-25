import requests
import json
from multiprocessing import Pool

NB_IDS_PER_REQUEST = 10
NB_PROCESSES = 16


def create_send_request(param : tuple):
    liste_ensembl, url, headers = param
    data = json.dumps({"ids": liste_ensembl})
    response = requests.post(url, headers=headers, data=data)
    if response.status_code == 200:
        print('Request successfull')
        return response.json()
    else:
        print(f"Failed to retrieve sequences: {response.status_code}, {response.text}")
        return []

        
# Fonction pour récupérer des séquences
def get_gene_sequences(ensembl_ids : list):
    url = "https://rest.ensembl.org/sequence/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    params = list()
    nb_id = len(ensembl_ids)
    gene_list = list()
    id = 0
    while id < nb_id :
        params.append((ensembl_ids[id:id+NB_IDS_PER_REQUEST], url, headers))
        id += NB_IDS_PER_REQUEST
    pool = Pool(processes = NB_PROCESSES)
    sequences = pool.map(create_send_request,params)
    for seq in sequences:
        gene_list.extend(seq)
    return gene_list



    
if __name__ == "__main__":
    ensembl_ids = [
    "ENSG00000139618", "ENSG00000141510", "ENSG00000157764", "ENSG00000198947", "ENSG00000198793",
    "ENSG00000198786", "ENSG00000198888", "ENSG00000198804", "ENSG00000198805", "ENSG00000198887"]

    # Récupérer les séquences en batchs (taille de batch = 10, 5 requêtes simultanées)
    import time
    debut = time.time()
    sequences = get_gene_sequences(ensembl_ids)

    # Affiche les premières séquences
    for seq in sequences:
        print(f"ID: {seq['id']}, Sequence: {seq['seq'][:100]}...\n")
    print(f"Temps écoulé : {time.time() - debut:.2f} secondes")
