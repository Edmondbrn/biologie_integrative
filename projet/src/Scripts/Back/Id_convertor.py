import mygene
import pandas as pd

def convert_refseq_to_ensembl(refseq_ids):
    """
    Convertit une liste d'identifiants RefSeq en identifiants Ensembl
    en utilisant la librairie MyGene.
    
    refseq_ids: liste de string (e.g. ["NM_001301415", "NR_002847", ...])
    Retour: liste de dictionnaires contenant les résultats.
    """
    mg = mygene.MyGeneInfo()
    # scopes='refseq' : on indique que l'on fournit des identifiants RefSeq
    # fields='ensembl.transcript' : on demande à récupérer les transcrits Ensembl
    # species='mouse' pour la souris (Mus musculus)
    results = mg.querymany(refseq_ids, 
                           scopes='refseq', 
                           fields='ensembl.transcript', 
                           species='mouse')
    return results

def add_ensembl_ids(file_path: str):
    df = pd.read_csv(file_path, sep='\t')
    df = df.drop(columns=["ensembl_transcript_id"])
    ncbi_ids = df["refseq"].tolist()
    conversion = convert_refseq_to_ensembl(ncbi_ids)
    
    # Créer un dictionnaire de conversion
    conversion_dict = {}
    for res in conversion:
        if 'ensembl' in res and 'transcript' in res['ensembl']:
            if not isinstance(res["ensembl"]["transcript"], list) :
                conversion_dict[res['query']] = res['ensembl']['transcript']
            else:
                conversion_dict[res['query']] = res['ensembl']['transcript'][0]
        else:
            conversion_dict[res['query']] = '\'Not Found\''
    
    # Ajouter une nouvelle colonne avec les identifiants Ensembl
    df['ensembl_id'] = df['refseq'].map(conversion_dict).fillna('Not Found')
    
    # Sauvegarder le fichier modifié
    df.to_csv(file_path +"_converted", sep='\t', index=False)

if __name__ == '__main__':
    file_path = "/home/edmond/Documents/GB5/biologie_integrative/projet/src/Ressources/data/FMRP_Binding_sites_mouse_Maurin_NAR_2014_merged.tsv"
    add_ensembl_ids(file_path)
