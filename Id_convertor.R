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
        write.table(merged_file, file = paste0("merged_", key, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
}
