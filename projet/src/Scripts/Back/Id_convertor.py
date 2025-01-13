import os
import pandas as pd
from pybiomart import Server

def check_available_filters():
    server = Server(host='http://www.ensembl.org')
    mart = server.marts['ENSEMBL_MART_ENSEMBL']
    dataset = mart.datasets['mmusculus_gene_ensembl']
    filters = dataset.filters
    print(filters)


# ------------------------------
# 1) Vérifier le type d'identifiants NCBI (NM, NR, XM, XR)
# ------------------------------
def check_data(identifier: str) -> str:
    """
    Retourne le type de RefSeq correspondant
    si l'identifiant contient 'NM', 'NR', 'XM' ou 'XR'.
    Lève une exception sinon.
    """
    if "NM" in identifier:
        return "refseq_mrna"
    elif "NR" in identifier:
        return "refseq_ncrna"
    elif "XM" in identifier:
        return "refseq_mrna_predicted"
    elif "XR" in identifier:
        return "refseq_ncrna_predicted"
    else:
        raise ValueError("Invalid identifier format. Please provide identifiers with NCBI prefixes (NM, NR, XM, XR).")

# ------------------------------
# 2) Conversion des identifiants NCBI vers Ensembl
# ------------------------------
def convert_ids(ncbi_ids, data_type):
    server = Server(host='http://www.ensembl.org')
    mart = server.marts['ENSEMBL_MART_ENSEMBL']
    dataset = mart.datasets['mmusculus_gene_ensembl']

    results = dataset.query(attributes=[data_type, 'ensembl_transcript_id'], filters={data_type: ncbi_ids})
    return results

# ------------------------------
# 3) Ajouter les identifiants Ensembl au fichier
# ------------------------------
def add_ensembl_ids(file_path: str):
    file_extension = os.path.splitext(file_path)[1].lower()
    
    # Lire le fichier en fonction de son extension
    if file_extension == '.csv':
        df = pd.read_csv(file_path)
    elif file_extension == '.tsv':
        df = pd.read_csv(file_path, sep='\t')
    elif file_extension in ['.xls', '.xlsx']:
        df = pd.read_excel(file_path)
    else:
        raise ValueError("Unsupported file format. Please provide a CSV, TSV, or Excel file.")
    
    # Identifier la colonne contenant les identifiants NCBI
    ncbi_column = None
    for column in df.columns:
        if any(prefix in str(df[column].iloc[0]) for prefix in ["NM", "NR", "XM", "XR"]):
            ncbi_column = column
            break
    
    if ncbi_column is None:
        raise ValueError("No NCBI identifiers found in the file.")
    
    # Obtenir tous les types de données NCBI présents dans la colonne
    data_types = df[ncbi_column].apply(check_data).unique()
    
    # Créer un dictionnaire de conversion pour chaque type de données NCBI
    conversion_dict = {}
    for data_type in data_types:
        ncbi_ids = df[df[ncbi_column].apply(check_data) == data_type][ncbi_column].astype(str).tolist()
        conversion_results = convert_ids(ncbi_ids, data_type)
        conversion_dict.update(dict(zip(conversion_results[data_type], conversion_results['ensembl_transcript_id'])))
    
    # Ajouter une nouvelle colonne avec les identifiants Ensembl
    df['ensembl_transcript_id'] = df[ncbi_column].map(conversion_dict).fillna('Not Found')
    
    # Sauvegarder le fichier modifié
    if file_extension == '.csv':
        df.to_csv(file_path, index=False)
    elif file_extension == '.tsv':
        df.to_csv(file_path, sep='\t', index=False)
    elif file_extension in ['.xls', '.xlsx']:
        df.to_excel(file_path, index=False)

