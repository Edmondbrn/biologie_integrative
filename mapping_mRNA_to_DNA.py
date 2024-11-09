import requests
from  pyensembl import EnsemblRelease
from multiprocessing import Pool
import time
import json

import aiohttp
import asyncio
from multiprocessing import Pool
from typing import List, Tuple, Dict

GENE_PER_PROCESS = 5
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
            sequence = json.loads(response.text) # convertit la str en dictionnaire
            liste_start, liste_end, chr, strand = [], [], "", ""
            for match in sequence["mappings"]: # browse all the matches
                start_gen = match["start"] # get the start and the end position
                end_gen = match["end"]
                if chr == "" and strand == "":
                    chr = match["seq_region_name"]
                    strand = match["strand"]

                new_start, new_end = ConvertCoordAssembly(start_gen, end_gen, specie, chr, strand)
                liste_start.append(new_start)
                liste_end.append(new_end)

            dict_sequences[transcrit_id] = {"start" : liste_start, "end" : liste_end}
            # TODO Voir pour diviser le cas où on a plusieurs fois le même transcrit
            # dict_sequences[transcrit_id] = response.text.strip()
        else:
            print(f"Failed to retrieve transcript ID: {response.status_code}")
            dict_sequences[transcrit_id] = {"start" : None, "end" : None}
    return dict_sequences

def ConvertCoordAssembly(start : int, end : int, specy : str, chromosome : int, strand : str, assembly1 : str  = "GRCm39", assembly2 : str = "GRCm38"):
    """
    Function to convert dna coordinate from one genome assmbly to another one
    """
    server = "https://rest.ensembl.org"
    ext = f"/map/{specy}/{assembly1}/{chromosome}:{start}..{end}:{strand}/{assembly2}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if r.status_code == 200:
        decoded = r.json()
        if decoded["mappings"] != []:
            update_start = decoded["mappings"][0]["mapped"]["start"]
            update_end = decoded["mappings"][0]["mapped"]["end"]
        else:
            update_start = None
            update_end = None
        return update_start, update_end

    else :
        print(f"Error cannot convert from {assembly1} to {assembly2}")
        return None, None
        

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
# import aiohttp
# import asyncio
# from typing import List, Tuple, Dict

# async def fetch(session, url):
#     async with session.get(url, headers={"Content-Type": "application/json"}) as response:
#         if response.headers['Content-Type'] == 'application/json':
#             return await response.json()
#         else:
#             text = await response.text()
#             print(f"Unexpected content type: {response.headers['Content-Type']}")
#             print(f"Response text: {text}")
#             return None

# async def request_gene_async(parameters: Tuple[List[str], List[int], List[int], str]):
#     trans_ids, start_list, end_list, specie = parameters
#     base_url = "https://rest.ensembl.org"
#     dict_sequences = {}

#     async with aiohttp.ClientSession() as session:
#         tasks = []
#         for transcrit_id, start, end in zip(trans_ids, start_list, end_list):
#             url = f"{base_url}/map/cdna/{transcrit_id}/{start}..{end}"
#             tasks.append(fetch(session, url))

#         responses = await asyncio.gather(*tasks)

#         for transcrit_id, response in zip(trans_ids, responses):
#             if response and "mappings" in response:
#                 liste_start, liste_end, chr, strand = [], [], "", ""
#                 for match in response["mappings"]:
#                     start_gen = match["start"]
#                     end_gen = match["end"]
#                     if chr == "" and strand == "":
#                         chr = match["seq_region_name"]
#                         strand = match["strand"]

#                     new_start, new_end = ConvertCoordAssembly(start_gen, end_gen, specie, chr, strand)
#                     liste_start.append(new_start)
#                     liste_end.append(new_end)

#                 dict_sequences[transcrit_id] = {"start": liste_start, "end": liste_end}
#             else:
#                 print(f"Failed to retrieve transcript ID: {transcrit_id}")
#                 dict_sequences[transcrit_id] = {"start": None, "end": None}

#     return dict_sequences

# def ConvertCoordAssembly(start: int, end: int, specy: str, chromosome: int, strand: str, assembly1: str = "GRCm39", assembly2: str = "GRCm38"):
#     server = "https://rest.ensembl.org"
#     ext = f"/map/{specy}/{assembly1}/{chromosome}:{start}..{end}:{strand}/{assembly2}?"
#     r = requests.get(server + ext, headers={"Content-Type": "application/json"})

#     if r.status_code == 200:
#         decoded = r.json()
#         if decoded["mappings"]:
#             update_start = decoded["mappings"][0]["mapped"]["start"]
#             update_end = decoded["mappings"][0]["mapped"]["end"]
#         else:
#             update_start = None
#             update_end = None
#         return update_start, update_end
#     else:
#         print(f"Error cannot convert from {assembly1} to {assembly2}")
#         return None, None

# def get_dna_sequences(transcript_ids_list: List[str], start_trans: List[int], end_trans: List[int], specy: str = "mouse"):
#     id = 0
#     nb_id = len(transcript_ids_list)
#     params = []
#     while id < nb_id:
#         params.append((transcript_ids_list[id:id + GENE_PER_PROCESS], start_trans[id:id + GENE_PER_PROCESS], end_trans[id:id + GENE_PER_PROCESS], specy))
#         id += GENE_PER_PROCESS

#     loop = asyncio.get_event_loop()
#     gene_sequences = loop.run_until_complete(asyncio.gather(*[request_gene_async(param) for param in params]))

#     final_dict = {}
#     for dict_genes in gene_sequences:
#         final_dict.update(dict_genes)
#     return final_dict

# # Exemple d'utilisation
# transcript_ids_list = [
#     "ENSMUST00000193812", "ENSMUST00000193813", "ENSMUST00000193814", "ENSMUST00000193815", "ENSMUST00000193816",
#     "ENSMUST00000193817", "ENSMUST00000193818", "ENSMUST00000193819", "ENSMUST00000193820", "ENSMUST00000193821",
#     "ENSMUST00000193822", "ENSMUST00000193823", "ENSMUST00000193824", "ENSMUST00000193825", "ENSMUST00000193826",
#     "ENSMUST00000193827", "ENSMUST00000193828", "ENSMUST00000193829", "ENSMUST00000193830", "ENSMUST00000193831"
# ]

# start_trans = [
#     1, 1, 1, 1, 1,
#     1, 1, 1, 1, 1,
#     1, 1, 1, 1, 1,
#     1, 1, 1, 1, 1
# ]

# end_trans = [
#     1000, 1100, 1200, 1300, 1400,
#     1500, 1600, 1700, 1800, 1900,
#     2000, 2100, 2200, 2300, 2400,
#     2500, 2600, 2700, 2800, 2900
# ]

# debut = time.time()
# # Appeler la fonction get_dna_sequences
# sequences = get_dna_sequences(transcript_ids_list, start_trans, end_trans, specy="mouse")

# # Afficher les séquences récupérées
# for transcript_id, coords in sequences.items():
#     print(f"Transcript ID: {transcript_id}, Start: {coords['start']}, End: {coords['end']}")
# print(f"Temps d'exécution: {time.time() - debut} secondes")
if __name__ == "__main__":
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
    debut = time.time()
    sequences = get_dna_sequences(transcript_ids_list, start_trans, end_trans, specy="mouse")

    # Afficher les séquences récupérées
    for transcript_id, dict_transcript in sequences.items():
        liste_start = dict_transcript["start"]
        liste_end = dict_transcript["end"]
        if liste_start != None or liste_end != None:

            print(f"Le transcrit {transcript_id} avec les coordonnées ARN {transcript_dict[transcript_id]} et génomique : {liste_start} - {liste_end} ")

        else:
            print(f"Transcript ID: {transcript_id} not found")
    print(f"Temps d'exécution: {time.time() - debut} secondes")