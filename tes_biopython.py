# import requests
# import json

# def get_gene_sequences_batch(ensembl_ids):
#     url = "https://rest.ensembl.org/sequence/id"
#     headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
#     # Prépare la requête en lot avec les IDs Ensembl
#     data = json.dumps({"ids": ensembl_ids})

#     # Requête POST
#     response = requests.post(url, headers=headers, data=data)

#     # Vérifie le statut de la requête
#     if response.status_code == 200:
#         sequences = response.json()
#         return sequences
#     else:
#         raise Exception(f"Failed to retrieve sequences: {response.status_code}, {response.text}")

# # Exemple avec 10 IDs Ensembl (remplacez par votre liste)
# ensembl_ids = [
#     "ENSG00000139618",  # BRCA2
#     "ENSG00000141510",  # TP53
#     "ENSG00000157764",  # BRAF
#     # Ajoutez jusqu'à 4000 IDs
# ]

# # Récupération des séquences pour tous les gènes
# sequences = get_gene_sequences_batch(ensembl_ids)

# # Affichage des résultats
# for seq in sequences:
#     print(f"ID: {seq['id']}\nSequence: {seq['seq'][:100]}...\n")

import requests
import json
from concurrent.futures import ThreadPoolExecutor

# Fonction pour récupérer des séquences en batch via POST
def get_gene_sequences_batch(ensembl_ids):
    url = "https://rest.ensembl.org/sequence/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = json.dumps({"ids": ensembl_ids})
    
    response = requests.post(url, headers=headers, data=data)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Failed to retrieve sequences: {response.status_code}, {response.text}")

# Fonction pour diviser la liste en sous-listes (batchs)
def batch_list(ensembl_ids, batch_size):
    for i in range(0, len(ensembl_ids), batch_size):
        yield ensembl_ids[i:i + batch_size]

# Liste de 100 IDs Ensembl (remplacez par vos propres IDs si nécessaire)
ensembl_ids = [
    "ENSG00000139618", "ENSG00000141510", "ENSG00000157764", "ENSG00000198947", "ENSG00000198793",
    "ENSG00000198786", "ENSG00000198888", "ENSG00000198804", "ENSG00000198805", "ENSG00000198887",
    "ENSG00000198886", "ENSG00000198899", "ENSG00000198900", "ENSG00000198794", "ENSG00000198795",
    "ENSG00000198792", "ENSG00000198789", "ENSG00000198809", "ENSG00000198810", "ENSG00000198811",
    "ENSG00000198812", "ENSG00000198813", "ENSG00000198814", "ENSG00000198815", "ENSG00000198816",
    "ENSG00000198817", "ENSG00000198818", "ENSG00000198819", "ENSG00000198820", "ENSG00000198821",
    "ENSG00000198822", "ENSG00000198823", "ENSG00000198824", "ENSG00000198825", "ENSG00000198826",
    "ENSG00000198827", "ENSG00000198828", "ENSG00000198829", "ENSG00000198830", "ENSG00000198831",
    "ENSG00000198832", "ENSG00000198833", "ENSG00000198834", "ENSG00000198835", "ENSG00000198836",
    "ENSG00000198837", "ENSG00000198838", "ENSG00000198839", "ENSG00000198840", "ENSG00000198841",
    "ENSG00000198842", "ENSG00000198843", "ENSG00000198844", "ENSG00000198845", "ENSG00000198846",
    "ENSG00000198847", "ENSG00000198848", "ENSG00000198849", "ENSG00000198850", "ENSG00000198851",
    "ENSG00000198852", "ENSG00000198853", "ENSG00000198854", "ENSG00000198855", "ENSG00000198856",
    "ENSG00000198857", "ENSG00000198858", "ENSG00000198859", "ENSG00000198860", "ENSG00000198861",
    "ENSG00000198862", "ENSG00000198863", "ENSG00000198864", "ENSG00000198865", "ENSG00000198866",
    "ENSG00000198867", "ENSG00000198868", "ENSG00000198869", "ENSG00000198870", "ENSG00000198871",
    "ENSG00000198872", "ENSG00000198873", "ENSG00000198874", "ENSG00000198875", "ENSG00000198876",
    "ENSG00000198877", "ENSG00000198878", "ENSG00000198879", "ENSG00000198880", "ENSG00000198881",
    "ENSG00000198882", "ENSG00000198883", "ENSG00000198884", "ENSG00000198885", "ENSG00000198886",
    "ENSG00000198887", "ENSG00000198888", "ENSG00000198889", "ENSG00000198890", "ENSG00000198891"
]

# Parallélisation des batchs
def fetch_all_sequences_in_batches(ensembl_ids, batch_size=10, max_workers=5):
    sequences = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        batch_generator = batch_list(ensembl_ids, batch_size)
        results = executor.map(get_gene_sequences_batch, batch_generator)
        
        for batch_result in results:
            sequences.extend(batch_result)
    
    return sequences

# Récupérer les séquences en batchs (taille de batch = 10, 5 requêtes simultanées)
sequences = fetch_all_sequences_in_batches(ensembl_ids, batch_size=10, max_workers=5)

# Affiche les premières séquences
for seq in sequences:
    print(f"ID: {seq['id']}, Sequence: {seq['seq'][:100]}...\n")
