import requests
import json
from multiprocessing import Pool
from pyensembl import EnsemblRelease

NB_IDS_PER_REQUEST = 50
NB_PROCESSES = 8
data = EnsemblRelease(species = 'mouse')


def create_send_request(param : tuple):
    liste_ensembl, url, headers = param
    data = json.dumps({"ids": liste_ensembl})
    response = session.post(url, headers=headers, data=data)
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
    nb_id = len(ensembl_ids)
    params = list()
    id = 0
    gene_list = list()
    while id < nb_id : # On découpe la liste des ensembl_ids en batchs de taille NB_IDS_PER_REQUEST
        params.append((ensembl_ids[id:id+NB_IDS_PER_REQUEST], url, headers))
        id += NB_IDS_PER_REQUEST

    pool = Pool(processes = NB_PROCESSES) # On crée un pool de processus pour la parallélisation
    sequences = pool.map(create_send_request,params) # on exécute les requêtes en parallèle
    for seq in sequences:
        gene_list.extend(seq)
    return gene_list
    


if __name__ == "__main__":
    session = requests.Session()

    ensembl_ids = [
    "ENSG00000139618", "ENSG00000141510", "ENSG00000157764", "ENSG00000198947", "ENSG00000198793",
    "ENSG00000198786", "ENSG00000198888", "ENSG00000198804", "ENSG00000198805", "ENSG00000198887"]

    # Récupérer les séquences en batchs (taille de batch = 10, 5 requêtes simultanées)
    import time
    debut = time.time()
    sequences = get_gene_sequences(ensembl_ids)

    # Affiche les premières séquences
    fh = open("test_ensembl_request.txt", "w")
    for seq in sequences:
        print(f"ID: {seq['id']}, Sequence: {seq['seq'][:100]}...\n")
        # save result in FASTA file avec retour à la ligne tous les 80 nulcetotides
        fh.write(f">{seq['id']}\n")
        for i in range(0, len(seq['seq']), 80):
            fh.write(seq['seq'][i:i+80] + "\n")
        fh.write("="*80 + "\n")

    fh.close()
    print(f"Temps écoulé : {time.time() - debut:.2f} secondes")



# import aiohttp
# import asyncio
# import json

# NB_IDS_PER_REQUEST = 10

# async def create_send_request(session, param):
#     liste_ensembl, url, headers = param
#     data = json.dumps({"ids": liste_ensembl})
#     async with session.post(url, headers=headers, data=data) as response:
#         if response.status == 200:
#             print('Request successful')
#             return await response.json()
#         else:
#             print(f"Failed to retrieve sequences: {response.status}, {await response.text()}")
#             return []

# async def get_gene_sequences(ensembl_ids):
#     url = "https://rest.ensembl.org/sequence/id"
#     headers = {"Content-Type": "application/json", "Accept": "application/json"}
#     params = [(ensembl_ids[i:i+NB_IDS_PER_REQUEST], url, headers) for i in range(0, len(ensembl_ids), NB_IDS_PER_REQUEST)]
    
#     async with aiohttp.ClientSession() as session:
#         tasks = [create_send_request(session, param) for param in params]
#         sequences = await asyncio.gather(*tasks)
    
#     gene_list = []
#     for seq in sequences:
#         gene_list.extend(seq)
#     return gene_list

# if __name__ == "__main__":
#     ensembl_ids = [
#         "ENSG00000139618", "ENSG00000141510", "ENSG00000157764", "ENSG00000198947", "ENSG00000198793",
#         "ENSG00000198786", "ENSG00000198888", "ENSG00000198804", "ENSG00000198805", "ENSG00000198887"
#     ]

#     # Récupérer les séquences en batchs (taille de batch = 10)
#     import time
#     debut = time.time()
#     sequences = asyncio.run(get_gene_sequences(ensembl_ids))

#     # Affiche les premières séquences
#     for seq in sequences:
#         print(f"ID: {seq['id']}, Sequence: {seq['seq'][:100]}...\n")
#     print(f"Temps écoulé : {time.time() - debut:.2f} secondes")