import requests
from  pyensembl import EnsemblRelease
from multiprocessing import Pool
import time
import json
from pandas import read_csv
import tqdm

GENE_PER_PROCESS = 20
NB_PROCESS = 16
REQUESTS_PER_SECOND = 15


 


def request_gene(parameters: tuple):
    trans_ids, start_list, end_list, specie = parameters
    # URL de base pour l'API REST d'Ensembl
    base_url = "https://rest.ensembl.org"
    dict_sequences = {}
    
    for transcrit_id, start, end in zip(trans_ids, start_list, end_list):#, total=len(trans_ids), desc="Processing transcripts"):
        while True:
            # Construire l'URL pour la requête
            url = f"{base_url}/map/cdna/{transcrit_id}/{start}..{end}"
            # Envoyer la requête GET
            response = requests.get(url, headers={"Content-Type": "application/json"})
            
            # Vérifier le statut de la réponse
            if response.status_code == 200:
                sequence = json.loads(response.text)  # convertit la str en dictionnaire
                liste_start, liste_end, chr, strand = [], [], "", ""
                for match in sequence["mappings"]:  # browse all the matches
                    start_gen = match["start"]  # get the start and the end position
                    end_gen = match["end"]
                    if chr == "" and strand == "":
                        chr = match["seq_region_name"]
                        strand = match["strand"]

                    new_start, new_end = ConvertCoordAssembly(start_gen, end_gen, specie, chr, strand)
                    liste_start.append(new_start)
                    liste_end.append(new_end)

                dict_sequences[transcrit_id] = {"start": liste_start, "end": liste_end}
                break
            elif response.status_code == 429:
                print("Rate limit exceeded. Waiting before retrying...")
                time.sleep(10)  # Attendre 10 secondes avant de réessayer
            else:
                print(f"Failed to retrieve transcript ID: {response.status_code}")
                dict_sequences[transcrit_id] = {"start": None, "end": None}
                break
    
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

if __name__ == "__main__":
    # Exemple d'utilisation
    df = read_csv("data_filtered2.tsv", sep = "\t", header = 0)
    df = df.loc[(df["start_ensembl"] != "unknown_id") | (df["end_ensembl"] != "Not found")]
    transcript_ids_list = df["ensembl_transcript_id"].tolist()

    start_trans = df["start_ensembl"].tolist()

    end_trans = df["end_ensembl"].tolist()
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