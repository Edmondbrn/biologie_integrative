import pandas as pd
import os

os.chdir(os.path.dirname(__file__))

def checkJCEC(a_file : str):
    return "JCEC" in a_file


def checkMATS(a_file : str):
    return "MATS" in a_file

def getRmatsFiles(rmats_path) -> list[str]:
    """
    Function to get only the 5 important files in rmats output
    """
    list_file = os.listdir(rmats_path)
    liste_to_keep = []
    for file in list_file:
        if checkJCEC(file) and  checkMATS(file):
            liste_to_keep.append(rmats_path + "/" + file)
    
    return liste_to_keep


def filterStrand(rmats_list):
    """
    Function to separate transcript on + and - strand
    """
    dict_file = {}
    for rmats_file in rmats_list: # browse the 5 important rmats files
        splicing = (rmats_file.split("/")[-1]).split(".")[0]
        df = pd.read_csv(rmats_file, sep ="\t")
        df_positive = df.loc[df["strand"] == "+"]
        df_negative = df.loc[df["strand"] == "-"]
        dict_file[f"{splicing}_+"] = df_positive
        dict_file[f"{splicing}_-"] = df_negative
    return dict_file


def filterA5SS(data : pd.DataFrame, splice_type : str) -> None:
    data = data[["GeneID", "chr", "strand", "shortEE", "longExonEnd", "flankingES"]]
    data.columns = ["GeneID", "chr", "strand", "shortSplice", "longSplice", "shareSplice"]
    data.to_csv(f"filteredRmats/{splice_type}_filtered.csv", sep = "\t")

def filterA3SS(data, splice_type) -> None:
    data = data[["GeneID", "chr", "strand", "shortES", "longExonStart_0base", "flankingEE"]]
    data.columns = ["GeneID", "chr", "strand", "shortSplice", "longSplice", "shareSplice"]
    data.to_csv(f"filteredRmats/{splice_type}_filtered.csv", sep = "\t")

def filterRI(data : pd.DataFrame, splice_type : str) -> None:
    data = data[["GeneID", "chr", "strand", "riExonStart_0base", "riExonEnd"]]
    data.columns = ["GeneID", "chr", "strand", "RiStart", "RiEnd"]
    data.to_csv(f"filteredRmats/{splice_type}_filtered.csv", sep = "\t")

def filterSE(data : pd.DataFrame, splice_type : str):
    data = data[["GeneID", "chr", "strand","exonStart_0base", "exonEnd", "upstreamEE", "downstreamES"]]
    data.columns = ["GeneID", "chr", "strand", "SEstart", "SeEnd", "upstreamEnd", "DownstreamStart"]
    data.to_csv(f"filteredRmats/{splice_type}_filtered.csv", sep = "\t")

def filterMXE(data : pd.DataFrame, splice_type : str) -> None:
    # Implémentez la logique de filtrage pour MXE+
    data = data[["GeneID", "chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamEE", "downstreamES"]]
    data.columns = ["GeneID", "chr", "strand", "1stExonStart", "1stExonEnd", "2ndExonStart", "2ndExonEnd", "upstreamEE", "downstreamES"]
    data.to_csv(f"filteredRmats/{splice_type}_filtered.csv", sep = "\t")
    pass



def chooseParsing(rmats_dict):
    """
    Function to choose what columns to keep for each type of splicing
    """
    dict_function = {
        "A5SS_+": filterA5SS,
        "A5SS_-": filterA3SS,
        "A3SS_+": filterA3SS,
        "A3SS_-": filterA5SS,
        "RI_+": filterRI,
        "RI_-": filterRI,
        "MXE_+": filterMXE,
        "MXE_-": filterMXE,
        "SE_+": filterSE,
        "SE_-": filterSE,
    }
    if not os.path.isdir("filteredRmats"):
        os.mkdir("filteredRmats")

    for element, data in rmats_dict.items():
        if element in dict_function:
            dict_function[element](data, element)  # Exécuter la fonction correspondante avec les données
            print(f"Filtered rmats file generated for {element}")



if __name__ == "__main__":
    liste_imp = getRmatsFiles("rmats_post")
    rmats_dico = filterStrand(liste_imp)
    chooseParsing(rmats_dico)