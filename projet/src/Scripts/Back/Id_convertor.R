# Charger le package biomaRt
library(biomaRt)

# Lire les identifiants NCBI à partir du fichier
check_data = function(file){
    if (grepl("NM", file, ignore.case = FALSE, perl = FALSE, fixed = TRUE))
        return("refseq_mrna")
    else if (grepl("NR", file, ignore.case = FALSE, perl = FALSE, fixed = TRUE))
        return("refseq_ncrna")
    else if (grepl("XM", file, ignore.case = FALSE, perl = FALSE, fixed = TRUE))
        return("refseq_mrna_predicted")
    else if (grepl("XR", file, ignore.case = FALSE, perl = FALSE, fixed = TRUE))
        return("refseq_ncrna_predicted")
    else
        stop("Invalid file format. Please provide a file with NCBI identifiers (NM, NR, XM, XR).")
}

convertor = function(file, data_type) {
    ncbi_ids <- read.table(file, header = FALSE, stringsAsFactors = FALSE)[, 1]
    # Connecter au serveur Biomart
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    # Effectuer la requête pour convertir les identifiants NCBI en identifiants Ensembl de transcrits
    results <- getBM(
        filters = data_type,
        attributes = c(data_type, "ensembl_transcript_id"),
        values = ncbi_ids,
        mart = ensembl
    )

    # Ajouter les noms des transcrits d'origine RefSeq dans le fichier de sortie
    return(results)
}

setwd("output")

results_list <- list(
    NM = list(),
    NR = list(),
    XM = list(),
    XR = list()
)

for (file in list.files()){
    print(file)
    data_type = check_data(file)
    result = convertor(file, data_type)
    if (data_type == "refseq_mrna")
        results_list$NM[[file]] = result
    else if (data_type == "refseq_ncrna")
        results_list$NR[[file]] = result
    else if (data_type == "refseq_mrna_predicted")
        results_list$XM[[file]] = result
    else if (data_type == "refseq_ncrna_predicted")
        results_list$XR[[file]] = result

}
# merging the NM files, XM files, NR files, XR files
for (key in names(results_list)){
    if (length(results_list[[key]]) > 0){
        merged_file = do.call(rbind, results_list[[key]])
        write.table(merged_file, file = paste0("merged_", key, ".tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
}
# suppression des fichiers intermédiaires
file.remove(list.files(pattern = ".txt"))


# ---------------------------------------------------------------------------------------------------------

# Code python à adapter pour produire le même effet

# ---------------------------------------------------------------------------------------------------------

# import os
# import glob
# from pybiomart import Server

# # ------------------------------
# # 1) Vérifier le type d'identifiants NCBI (NM, NR, XM, XR)
# # ------------------------------
# def check_data(file_name: str) -> str:
#     """
#     Retourne le type de RefSeq correspondant
#     si le nom de fichier contient 'NM', 'NR', 'XM' ou 'XR'.
#     Lève une exception sinon.
#     """
#     if "NM" in file_name:
#         return "refseq_mrna"
#     elif "NR" in file_name:
#         return "refseq_ncrna"
#     elif "XM" in file_name:
#         return "refseq_mrna_predicted"
#     elif "XR" in file_name:
#         return "refseq_ncrna_predicted"
#     else:
#         raise ValueError("Invalid file format. Please provide a file with NCBI identifiers (NM, NR, XM, XR).")

# # ------------------------------
# # 2) Conversion des identifiants NCBI vers Ensembl
# # ------------------------------
# def convertor(file_path: str, data_type: str, dataset) -> 'pd.DataFrame':
#     """
#     Lit la liste d'identifiants NCBI depuis le fichier
#     et utilise pybiomart pour récupérer les ID transcrits Ensembl correspondants.
#     """
#     # Lecture des identifiants (une ligne = un ID)
#     with open(file_path, 'r') as f:
#         ncbi_ids = [line.strip() for line in f if line.strip()]

#     # Requête via pybiomart
#     # On utilise le paramètre 'filters' de la forme { data_type: [ ... ] }
#     # et on demande les colonnes d'intérêt via le paramètre 'attributes'.
#     result_df = dataset.query(
#         attributes=[data_type, 'ensembl_transcript_id'],
#         filters={data_type: ncbi_ids}
#     )
#     return result_df

# # ------------------------------
# # 3) Script principal
# # ------------------------------
# if __name__ == "__main__":
    
#     # Récupération du serveur ensembl
#     server = Server(host='http://www.ensembl.org')
#     # On choisit la base "ENSEMBL_MART_ENSEMBL" et l'espèce Mus musculus
#     dataset = server['ENSEMBL_MART_ENSEMBL']['mmusculus_gene_ensembl']

#     # Passage au répertoire "output"
#     os.chdir("output")

#     # Dictionnaire pour stocker les dataframes de résultats
#     results_list = {
#         'NM': [],
#         'NR': [],
#         'XM': [],
#         'XR': []
#     }

#     # Pour chaque fichier présent dans le dossier courant
#     for file_name in glob.glob("*"):
#         # Éviter de traiter notre script Python, .pyc ou autres fichiers indésirables
#         if not os.path.isfile(file_name):
#             continue
        
#         print(f"Traitement du fichier : {file_name}")
        
#         # Détermine le type de données (ex: 'refseq_mrna')
#         data_type = check_data(file_name)

#         # Conversion
#         try:
#             result = convertor(file_name, data_type, dataset)
#         except Exception as e:
#             print(f"Erreur lors de la conversion pour {file_name} : {e}")
#             continue
        
#         # On range le DataFrame dans la liste correspondant à son type
#         if data_type == "refseq_mrna":
#             results_list['NM'].append(result)
#         elif data_type == "refseq_ncrna":
#             results_list['NR'].append(result)
#         elif data_type == "refseq_mrna_predicted":
#             results_list['XM'].append(result)
#         elif data_type == "refseq_ncrna_predicted":
#             results_list['XR'].append(result)

#     # ------------------------------
#     # 4) Fusion et écriture des résultats
#     # ------------------------------
#     import pandas as pd

#     for key, list_of_dfs in results_list.items():
#         if list_of_dfs:
#             # Fusion des dataframes pour le type `key`
#             merged_df = pd.concat(list_of_dfs, ignore_index=True)
#             # Écriture dans un TSV
#             output_file = f"merged_{key}.tsv"
#             merged_df.to_csv(output_file, sep="\t", index=False)
#             print(f"Fichier généré : {output_file}")

#     # ------------------------------
#     # 5) Suppression des fichiers intermédiaires
#     #    (ici, tous les .txt du répertoire)
#     # ------------------------------
#     for txt_file in glob.glob("*.txt"):
#         try:
#             os.remove(txt_file)
#             print(f"Fichier supprimé : {txt_file}")
#         except Exception as e:
#             print(f"Impossible de supprimer {txt_file} : {e}")

