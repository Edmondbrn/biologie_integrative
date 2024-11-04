import requests
from multiprocessing import Pool
import json


GENE_PER_PROCESS = 10
NB_PROCESS = 8
REQUESTS_PER_SECOND = 15

def request_gene(parameters: tuple):
    trans_ids, start_list, end_list, specie = parameters
    base_url = "https://rest.ensembl.org"
    dict_sequences = {}
    for transcrit_id, start, end in zip(trans_ids, start_list, end_list):
        url = f"{base_url}/map/cdna/{transcrit_id}/{start}..{end}"
        response = requests.get(url, headers={"Content-Type": "application/json"})
        if response.status_code == 200:
            dict_sequences[transcrit_id] = response.text.strip()
        else:
            print(f"Failed to retrieve transcript ID: {response.status_code}")
            dict_sequences[transcrit_id] = None
    return dict_sequences

def get_dna_sequences(transcript_ids_list: list[str], start_trans: list[int], end_trans: list[int], specy: str = "mouse"):
    id = 0
    nb_id = len(transcript_ids_list)
    params = []
    while id < nb_id:
        params.append((transcript_ids_list[id:id + GENE_PER_PROCESS], start_trans[id:id + GENE_PER_PROCESS], end_trans[id:id + GENE_PER_PROCESS], specy))
        id += GENE_PER_PROCESS

    with Pool(processes=NB_PROCESS) as pool:
        gene_sequences = pool.map(request_gene, params)
    final_dict = {}
    for dict_genes in gene_sequences:
        final_dict = {**final_dict, **dict_genes}
    return final_dict

if __name__ == "__main__":
    transcript_ids_list = [
        "ENSMUST00000193812", "ENSMUST00000193813", "ENSMUST00000193814", "ENSMUST00000193815", "ENSMUST00000193816",
        "ENSMUST00000193817", "ENSMUST00000193818", "ENSMUST00000193819", "ENSMUST00000193820", "ENSMUST00000193821",
        "ENSMUST00000193822", "ENSMUST00000193823", "ENSMUST00000193824", "ENSMUST00000193825", "ENSMUST00000193826",
        "ENSMUST00000193827", "ENSMUST00000193828", "ENSMUST00000193829", "ENSMUST00000193830", "ENSMUST00000193831"
    ]

    start_trans = [1] * len(transcript_ids_list)
    end_trans = [1000 + i * 100 for i in range(len(transcript_ids_list))]
    transcript_dict = {transcript_id: (start, end) for transcript_id, start, end in zip(transcript_ids_list, start_trans, end_trans)}

    sequences = get_dna_sequences(transcript_ids_list, start_trans, end_trans, specy="mouse")
    for transcript_id, sequence in sequences.items():
        if sequence is not None:
            sequence = json.loads(sequence)
            start_gen = sequence["mappings"][0]["start"]
            start_arn = transcript_dict[transcript_id][0]
            end_gen = sequence["mappings"][0]["end"]
            end_arn = transcript_dict[transcript_id][1]
            print(f"Transcript ID: {transcript_id} position {start_arn}:{end_arn}, Genomic position: {start_gen}:{end_gen}")
        else:
            print(f"Transcript ID: {transcript_id} not found")
