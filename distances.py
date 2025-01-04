import os
import pandas as pd
import numpy as np
import ast
import pyensembl as pb
from numba import njit
from multiprocessing import Pool, cpu_count
import pyensembl as pb
from distances_utils import convert_dna_to_rna

os.chdir(os.path.dirname(__file__))

class Distances():


    ERROR_DICT = {4 : "Error while converting dna to rna",
                    1 : "Not on the same transcript",
                    2 : "Second coordinate not in the transcript",
                    3 : "Error while searching for the second coordinates",
                    0 : "No problemo"}

    """
    This class is for computing distances between splicing sites and protein fixation site on mRNA
    """


    def __init__(this, ensembl_release: int = 102, specy: str = "mus_musculus"):
        this.__data_prot: pd.DataFrame = pd.DataFrame()
        this.__data_splicing: dict[str, pd.DataFrame] = dict()
        this.bdd = pb.EnsemblRelease(release = ensembl_release, species = specy)  # Déclaration de la variable de classe

    def LoadProgress(this, total: int) -> None:
        """
        Méthod to load  a progress bar.
        """
        this.progressbar = Progress(total)
        this.progressbar.show()

    def _LoadDataProt(this, file_name : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data of the protein
        file_name : name of the file of the protein
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        this.__data_prot = pd.read_csv(file_name, sep = sep)
        this.__data_prot["start_genomic"] = this.__data_prot["start_genomic"].apply(ast.literal_eval)
        this.__data_prot["end_genomic"] = this.__data_prot["end_genomic"].apply(ast.literal_eval)
        this.__data_prot = this.__FilterDataProt(this.__data_prot)

    def __FilterDataProt(this, df_prot : pd.DataFrame) -> pd.DataFrame:
        # enleve les lignes avec des str sur la colonne start_ensembl
        df_prot = df_prot.loc[df_prot["start_ensembl"].apply(lambda x: x.isnumeric())]
        return df_prot

    
    def __LoadDataSplicing(this, path : str, sep : str = "\t") -> pd.DataFrame:
        """
        Method to load the data from alternative splicing
        path : path to the folder containing the files
        sep : separator in the csv file [e.g : \t, `,` , ;]
        """
        for file in os.listdir(path):
            this.__data_splicing[file.split(".")[0]] = pd.read_csv(f"{path}/{file}", sep = sep)
    
    def __IsDataFrameNull(this, df : pd.DataFrame) -> bool:
        """
        Method to check if a dataframe is null
        """
        return df.shape[0] == 0


    def __CreateDistanceFile(this, df : pd.DataFrame, splice_type : str, output_dir : str, file_basename : str = "distances", type : str = "DNA") -> None:
        """
        Method to create a csv file containing the distances
        """
        if os.path.exists(output_dir) == False:
            os.mkdir(output_dir)
        df.to_csv(f"{output_dir}/{file_basename}_{type}_{splice_type}.csv", index = False, sep = "\t")

    @staticmethod
    @njit(fastmath = True)
    def __compute_distances(start_genomic_first : int, 
                            end_genomic_last : int, 
                            short_splice : int, 
                            share_splice : int, 
                            long_splice : int, 
                            exon_pos_list : list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif. 
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour A5SS et A3SS
        """
        # Initialisation des tableaux numpy pour stocker les résultats
        distances = np.zeros(12, dtype=np.int64)  # 12 distances
        flag = np.zeros(12, dtype=np.bool_)
        err_message = np.zeros(12, dtype=np.int64)

        distances[0] = start_genomic_first - short_splice
        distances[1] = start_genomic_first - share_splice
        distances[2] = end_genomic_last - short_splice
        distances[3] = end_genomic_last - share_splice
        distances[4] = start_genomic_first - long_splice
        distances[5] = end_genomic_last - long_splice
        distances[6], flag[6], err_message[6] = convert_dna_to_rna(start_genomic_first, short_splice, distances[0], exon_pos_list)
        distances[7], flag[7], err_message[7]  = convert_dna_to_rna(start_genomic_first, share_splice, distances[1], exon_pos_list)
        distances[8], flag[8], err_message[8]  = convert_dna_to_rna(end_genomic_last, short_splice, distances[2], exon_pos_list)
        distances[9], flag[9], err_message[9]  = convert_dna_to_rna(end_genomic_last, share_splice, distances[3], exon_pos_list)
        distances[10], flag[10], err_message[10]  = convert_dna_to_rna(start_genomic_first, long_splice, distances[4], exon_pos_list)
        distances[11], flag[11], err_message[11]  = convert_dna_to_rna(end_genomic_last, long_splice, distances[5], exon_pos_list)
        return distances, flag, err_message


    @staticmethod
    @njit(fastmath = True)
    def __compute_distances_RI(start_genomic_first : int, 
                               end_genomic_last : int, 
                               RiStart : int, 
                               RiEnd : int, 
                               exon_pos_list : list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif. 
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour la rétention d'intron
        """
        # Initialisation des tableaux numpy pour stocker les résultats
        distances = np.zeros(8, dtype=np.int64)  # 4 distances
        flag = np.zeros(8, dtype=np.bool_)
        err_message = np.zeros(8, dtype=np.int64)

        distances[0] = start_genomic_first - RiStart
        distances[1] = start_genomic_first - RiEnd
        distances[2] = end_genomic_last - RiStart
        distances[3] = end_genomic_last - RiEnd
        distances[4], flag[4], err_message[4] = convert_dna_to_rna(start_genomic_first, RiStart, distances[0], exon_pos_list)
        distances[5], flag[5], err_message[5] = convert_dna_to_rna(start_genomic_first, RiEnd, distances[1], exon_pos_list)
        distances[6], flag[6], err_message[6] = convert_dna_to_rna(end_genomic_last, RiStart, distances[2], exon_pos_list)
        distances[7], flag[7], err_message[7] = convert_dna_to_rna(end_genomic_last, RiEnd, distances[3], exon_pos_list)
        return distances, flag, err_message
    
    @staticmethod
    @njit(fastmath = True)
    def __compute_distances_SE(start_genomic_first : int, 
                               end_genomic_last : int,
                               upstreamEnd : int, 
                               DownstreamStart : int,  
                               exon_pos_list : list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif.
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour l'exon sauté (skipping exon)
        """
        # Initialisation des tableaux numpy pour stocker les résultats
        distances = np.zeros(8, dtype=np.int64)
        flag = np.zeros(8, dtype=np.bool_)
        err_message = np.zeros(8, dtype=np.int64)

        distances[0] = start_genomic_first - upstreamEnd
        distances[1] = start_genomic_first - DownstreamStart
        distances[2] = end_genomic_last - upstreamEnd
        distances[3] = end_genomic_last - DownstreamStart
        distances[4], flag[4], err_message[4] = convert_dna_to_rna(start_genomic_first, upstreamEnd, distances[0], exon_pos_list)
        distances[5], flag[5], err_message[5] = convert_dna_to_rna(start_genomic_first, DownstreamStart, distances[1], exon_pos_list)
        distances[6], flag[6], err_message[6] = convert_dna_to_rna(end_genomic_last, upstreamEnd, distances[2], exon_pos_list)
        distances[7], flag[7], err_message[7] = convert_dna_to_rna(end_genomic_last, DownstreamStart, distances[3], exon_pos_list)
        return distances, flag, err_message

    
    @staticmethod
    @njit(fastmath = True)
    def __compute_distances_MXE(start_genomic_first : int , 
                                end_genomic_last : int , 
                                FirstExonStart : int, 
                                FirstExonEnd : int, 
                                SecondExonStart: int, 
                                SecondExonEnd : int, 
                                upstreamEE : int, 
                                downstreamES : int, 
                                exon_pos_list : list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Calcule les distances pour un couple site de fixation de la protéine / site d'épissage alternatif.
        les arguments sont les différents coordonnées génomiques des acteurs
        Spécifique pour mutually exclusive exon
        """
        # Initialisation des tableaux numpy pour stocker les résultats
        distances = np.zeros(24, dtype=np.int64)  # 16 distances
        flag = np.zeros(24, dtype=np.bool_)
        err_message = np.zeros(24, dtype=np.int64)
        # Calcul des distances ADN
        distances[0] = start_genomic_first - upstreamEE
        distances[1] = start_genomic_first - FirstExonStart
        distances[2] = start_genomic_first - FirstExonEnd
        distances[3] = start_genomic_first - downstreamES
        distances[4] = end_genomic_last - upstreamEE
        distances[5] = end_genomic_last - FirstExonStart
        distances[6] = end_genomic_last - FirstExonEnd
        distances[7] = end_genomic_last - downstreamES
        distances[8] = start_genomic_first - SecondExonStart
        distances[9] = start_genomic_first - SecondExonEnd
        distances[10] = end_genomic_last - SecondExonStart
        distances[11] = end_genomic_last - SecondExonEnd
        # Conversion des distances ADN en distances ARN
        distances[12], flag[12], err_message[12] = convert_dna_to_rna(start_genomic_first, upstreamEE, distances[0], exon_pos_list)
        distances[13], flag[13], err_message[13] = convert_dna_to_rna(start_genomic_first, FirstExonStart, distances[1], exon_pos_list)
        distances[14], flag[14], err_message[14] = convert_dna_to_rna(start_genomic_first, FirstExonEnd, distances[2], exon_pos_list)
        distances[15], flag[15], err_message[15] = convert_dna_to_rna(start_genomic_first, downstreamES, distances[3], exon_pos_list)
        distances[16], flag[16], err_message[16] = convert_dna_to_rna(end_genomic_last, upstreamEE, distances[4], exon_pos_list)
        distances[17], flag[17], err_message[17] = convert_dna_to_rna(end_genomic_last, FirstExonStart, distances[5], exon_pos_list)
        distances[18], flag[18], err_message[18] = convert_dna_to_rna(end_genomic_last, FirstExonEnd, distances[6], exon_pos_list)
        distances[19], flag[19], err_message[19] = convert_dna_to_rna(end_genomic_last, downstreamES, distances[7], exon_pos_list)
        distances[20], flag[20], err_message[20] = convert_dna_to_rna(start_genomic_first, SecondExonStart, distances[8], exon_pos_list)
        distances[21], flag[21], err_message[21] = convert_dna_to_rna(start_genomic_first, SecondExonEnd, distances[9], exon_pos_list)
        distances[22] , flag[22], err_message[22] = convert_dna_to_rna(end_genomic_last, SecondExonStart, distances[10], exon_pos_list)
        distances[23] , flag[23], err_message[23] = convert_dna_to_rna(end_genomic_last, SecondExonEnd, distances[11], exon_pos_list)
        return distances, flag, err_message

    def __fill_rna_row(this, rna_indices : dict, dist_array : int, flag : bool, err_message :int , 
                       transcript_id : int, protein_sequence : str):
        """
        Méthode pour remplir le slignes des différents tableaux pour les distances ARN
        """
        row_rna = {}
        err_message = err_message[len(err_message)//2:]  # on ne garde que la 2e moitié qui correspond aux distances ARN
        for index, (i, col_name) in enumerate(rna_indices.items()):
            if err_message[index] != 0:
                # Remplacer la valeur par le message d’erreur ou un code
                row_rna[col_name] = f"ERROR_{Distances.ERROR_DICT[err_message[index]]}"
            elif flag[i]:
                # Ajouter un astérisque
                row_rna[col_name] = f"{dist_array[i]}*"
            else:
                # Valeur normale
                row_rna[col_name] = dist_array[i]
        row_rna["transcript_ID"]= transcript_id
        row_rna["prot_seq"]= protein_sequence
        return row_rna


    def distanceA5SS(this, splice_type: str = "A5SS_+",  outputdir : str = "distances/", basename : str = "distances_MXE") -> None:
        """
        Method which compute and format all the distances for the A5SS splicing
        It will create a csv file with the distances
        splice_type : type of splicing with the strand
        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]]
        results_dna = []  # On stockera ici toutes les lignes calculées
        results_ran = []
        for i in range(this.__data_prot.shape[0]): # browse all the protein fixation sites
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]] # Select the splicing sites of the same gene to avoid heavy computational loss
            if this.__IsDataFrameNull(splicing_same_gene):
                continue

            # Récupération des valeurs start/end (première et dernière coordonnée de la protéine)
            start_first = this.__data_prot.loc[i, "start_genomic"][0]
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]

            rna_indices = { # name of the columns in the final DataFrame
                6:  "prot_start_short_splice_start",
                7:  "prot_start_downstream_start",
                8:  "prot_end_short_splice_start",
                9:  "prot_end_downstream_start",
                10: "prot_start_long_splice_start",
                11: "prot_end_long_splice_start"
            }
            
            for j in range(splicing_same_gene.shape[0]): # browse all the splicing sites of the same gene
                short_sp = splicing_same_gene.iloc[j]["shortSplice"]
                long_sp = splicing_same_gene.iloc[j]["longSplice"]
                share_sp = splicing_same_gene.iloc[j]["shareSplice"]
                
                transcript_id = this.__data_prot.loc[i, "ensembl_id"]
                transcript : pb.Transcript = this.bdd.transcript_by_id(transcript_id) # get the transcript by the id to get exon informatinon"
                exon_pos_list = transcript.exon_intervals
                # Calcul des distances
                dist_array, flag, err_message = Distances.__compute_distances(start_first, end_last, short_sp, share_sp, long_sp, exon_pos_list)
                
                # Construire un dict pour tout remettre dans un DataFrame final
                row_dna = {
                    "prot_start_short_splice_start": dist_array[0], "prot_start_downstream_start":   dist_array[1], "prot_end_short_splice_start":   dist_array[2],
                    "prot_end_downstream_start":     dist_array[3], "prot_start_long_splice_start":  dist_array[4], "prot_end_long_splice_start":    dist_array[5],
                    "transcript_ID": transcript_id, "prot_seq": this.__data_prot.loc[i, "seq"]
                }

                row_rna = this.__fill_rna_row(rna_indices, dist_array, flag, err_message, transcript_id, this.__data_prot.loc[i, "seq"])

                results_dna.append(row_dna)
                results_ran.append(row_rna)
        # Convertir la liste de dicts en DataFrame
        df_dna = pd.DataFrame(results_dna)
        df_rna = pd.DataFrame(results_ran)
        this.__CreateDistanceFile(df_dna, splice_type,outputdir, basename)
        this.__CreateDistanceFile(df_rna, splice_type,outputdir, basename ,"RNA")




    def distanceA3SS(this, splice_type : str = "A3SS_+",  outputdir : str = "distances/", basename : str = "distances_MXE") :
        """
        Method which compute and format all the distances for the A3SS splicing
        It will create a csv file with the distances
        splice_type : type of splicing with the strand
        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "shortSplice", "longSplice", "shareSplice"]]
        results_dna = []  # On stockera ici toutes les lignes calculées
        results_rna = []
        for i in range(this.__data_prot.shape[0]):
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            # Récupération des valeurs start/end (première et dernière coordonnée)
            start_first = this.__data_prot.loc[i, "start_genomic"][0]
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]

            rna_indices = {
                6:  "prot_start_short_splice_start",
                7:  "prot_start_upstream_end",
                8:  "prot_end_short_splice_start",
                9:  "prot_end_upstream_end",
                10: "prot_start_long_splice_start",
                11: "prot_end_long_splice_start"
            }
            
            for j in range(splicing_same_gene.shape[0]):
                short_sp = splicing_same_gene.iloc[j]["shortSplice"]
                long_sp = splicing_same_gene.iloc[j]["longSplice"]
                share_sp = splicing_same_gene.iloc[j]["shareSplice"]
                transcript_id = this.__data_prot.loc[i, "ensembl_id"]
                transcript : pb.Transcript = this.bdd.transcript_by_id(transcript_id)
                exon_pos_list = transcript.exon_intervals

                dist_array, flag, err_message = Distances.__compute_distances(start_first, end_last, short_sp, share_sp, long_sp, exon_pos_list)
                # Construire un dict pour tout remettre dans un DataFrame final
                row_dna = {
                    "prot_start_short_splice_start": dist_array[0], "prot_start_upstream_end":   dist_array[1], "prot_end_short_splice_start":   dist_array[2],
                    "prot_end_upstream_end":     dist_array[3], "prot_start_long_splice_start":  dist_array[4], "prot_end_long_splice_start":    dist_array[5],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"], "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                row_rna = this.__fill_rna_row(rna_indices, dist_array, flag, err_message, transcript_id, this.__data_prot.loc[i, "seq"])
                results_dna.append(row_dna)
                results_rna.append(row_rna)
        # Convertir la liste de dicts en DataFrame
        df_dna = pd.DataFrame(results_dna)
        df_rna = pd.DataFrame(results_rna)
        this.__CreateDistanceFile(df_dna, splice_type,outputdir, basename)
        this.__CreateDistanceFile(df_rna, splice_type,outputdir, basename ,"RNA")


    def distanceRI(this, splice_type: str = "RI_+",  outputdir : str = "distances/", basename : str = "distances_MXE") -> None:
        """
        Method which compute and format all the distances for the RI splicing
        It will create a csv file with the distances
        splice_type : type of splicing with the strand
        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "RiStart", "RiEnd"]]  # récupère les coordonnées des sites de splicing
        results_dna = []
        results_rna = list()
        # Browse all the protein fixation sites
        for i in range(this.__data_prot.shape[0]):
            # Select the splicing sites of the same gene to avoid heavy computational loss
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            start_first = this.__data_prot.loc[i, "start_genomic"][0]
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]

            rna_indices = {
                4:  "prot_start_RiStart",
                5:  "prot_start_RiEnd",
                6:  "prot_end_RiStart",
                7:  "prot_end_RiEnd"
            }
            for j in range(splicing_same_gene.shape[0]):
                RiStart = splicing_same_gene.iloc[j]["RiStart"]
                RiEnd = splicing_same_gene.iloc[j]["RiEnd"]
                transcript_id = this.__data_prot.loc[i, "ensembl_id"]
                transcript = this.bdd.transcript_by_id(transcript_id)
                exon_pos_list = transcript.exon_intervals
                dist_array, flag, err_message = Distances.__compute_distances_RI(start_first, end_last, RiStart, RiEnd, exon_pos_list)
                row_dna = {
                    "prot_start_RiStart": dist_array[0], "prot_start_RiEnd": dist_array[1], "prot_end_RiStart": dist_array[2],
                    "prot_end_RiEnd": dist_array[3], "transcript_ID": this.__data_prot.loc[i, "ensembl_id"], "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                row_rna = this.__fill_rna_row(rna_indices, dist_array, flag, 
                                              err_message, transcript_id, this.__data_prot.loc[i, "seq"])
                results_dna.append(row_dna)
                results_rna.append(row_rna)
        df_dna = pd.DataFrame(results_dna)
        df_rna = pd.DataFrame(results_rna)
        this.__CreateDistanceFile(df_dna, splice_type,outputdir, basename)
        this.__CreateDistanceFile(df_rna, splice_type,outputdir, basename ,"RNA")
        
    def distanceSE(this, splice_type: str = "SE",  outputdir : str = "distances/", basename : str = "distances_MXE") -> None:
        """
        Method which compute and format all the distances for the SE splicing
        It will create a csv file with the distances
        splice_type : type of splicing with the strand
        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "upstreamEnd", "DownstreamStart"]]
        results_dna = []
        results_rna = []
        # Browse all the protein fixation sites
        for i in range(this.__data_prot.shape[0]):
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            start_first = this.__data_prot.loc[i, "start_genomic"][0] # extrait les coordonnées de fixation de la protéine
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]
            rna_indices = {
                4:  "prot_start_upstreamEnd",
                5:  "prot_start_DownstreamStart",
                6:  "prot_end_upstreamEnd",
                7:  "prot_end_DownstreamStart"
            }
            for j in range(splicing_same_gene.shape[0]):
                upstreamEnd = splicing_same_gene.iloc[j]["upstreamEnd"]
                DownstreamStart = splicing_same_gene.iloc[j]["DownstreamStart"]
                transcript_id = this.__data_prot.loc[i, "ensembl_id"]
                transcript : pb.Transcript = this.bdd.transcript_by_id(transcript_id)
                exon_pos_list = transcript.exon_intervals
                dist_array, flag, err_message = Distances.__compute_distances_SE(start_first, end_last, upstreamEnd, DownstreamStart, exon_pos_list)
                row_dna = {
                    "prot_start_upstreamEnd": dist_array[0], "prot_start_DownstreamStart": dist_array[1], "prot_end_upstreamEnd": dist_array[2],
                    "prot_end_DownstreamStart": dist_array[3], "transcript_ID": this.__data_prot.loc[i, "ensembl_id"], "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                row_rna = this.__fill_rna_row(rna_indices, dist_array, flag, err_message, transcript_id, this.__data_prot.loc[i, "seq"])
                results_dna.append(row_dna)
                results_rna.append(row_rna)
        df_dna = pd.DataFrame(results_dna)
        df_rna = pd.DataFrame(results_rna)
        this.__CreateDistanceFile(df_dna, splice_type,outputdir, basename)
        this.__CreateDistanceFile(df_rna, splice_type,outputdir, basename ,"RNA")

    def distanceMXE(this, splice_type: str = "MSE", outputdir : str = "distances/", basename : str = "distances_MXE") -> None:
        """
        Fonction pour obtenir les coordonnées des sites de splicing et les distances avec les sites de fixation des protéines
        splice_type : type de splicing
        """
        splicing_coordinates: pd.DataFrame = this.__data_splicing[splice_type][["GeneID", "chr", "strand", "1stExonStart", "1stExonEnd", "2ndExonStart", "2ndExonEnd", "upstreamEE", "downstreamES"]]
        results_dna = []
        results_rna = []
        # Browse all the protein fixation sites
        for i in range(this.__data_prot.shape[0]):
            splicing_same_gene = splicing_coordinates.loc[splicing_coordinates["GeneID"] == this.__data_prot.loc[i, "GeneID"]]
            if this.__IsDataFrameNull(splicing_same_gene):
                continue
            start_first = this.__data_prot.loc[i, "start_genomic"][0] # extrait les coordonnées de fixation de la protéine
            end_last = this.__data_prot.loc[i, "end_genomic"][-1]
            rna_indices = {
                12:  "prot_start_upstreamEE",
                13:  "prot_start_FirstExonStart",
                14:  "prot_start_FirstExonEnd",
                15:  "prot_start_downstreamES",
                16: "prot_end_upstreamEE",
                17: "prot_end_FirstExonStart",
                18 : "prot_end_FirstExonEnd",
                19 : "prot_end_downstreamES",
                20 : "prot_start_SecondExonStart",
                21 : "prot_start_SecondExonEnd",
                22 : "prot_end_SecondExonStart",
                23 : "prot_end_SecondExonEnd"
            }
            for j in range(splicing_same_gene.shape[0]):
                FirstExonStart = splicing_same_gene.iloc[j]["1stExonStart"]
                FirstExonEnd = splicing_same_gene.iloc[j]["1stExonEnd"]
                SecondExonStart = splicing_same_gene.iloc[j]["2ndExonStart"]
                SecondExonEnd = splicing_same_gene.iloc[j]["2ndExonEnd"]
                upstreamEE = splicing_same_gene.iloc[j]["upstreamEE"]
                downstreamES = splicing_same_gene.iloc[j]["downstreamES"]
                transcript_id = this.__data_prot.loc[i, "ensembl_id"]
                transcript = this.bdd.transcript_by_id(transcript_id)
                exon_pos_list = transcript.exon_intervals
                dist_array, flag, err_message = Distances.__compute_distances_MXE(start_first, end_last, FirstExonStart, FirstExonEnd, SecondExonStart, SecondExonEnd, upstreamEE, downstreamES, exon_pos_list)
                row_dna = {
                    "prot_start_upstreamEE": dist_array[0], "prot_start_FirstExonStart": dist_array[1], "prot_start_FirstExonEnd": dist_array[2],
                    "prot_start_downstreamES": dist_array[3], "prot_end_upstreamEE": dist_array[4], "prot_end_FirstExonStart": dist_array[5],
                    "prot_end_FirstExonEnd": dist_array[6], "prot_end_downstreamES": dist_array[7], "prot_start_SecondExonStart": dist_array[8],
                    "prot_start_SecondExonEnd": dist_array[9], "prot_end_SecondExonStart": dist_array[10], "prot_end_SecondExonEnd": dist_array[11],
                    "transcript_ID": this.__data_prot.loc[i, "ensembl_id"], "prot_seq": this.__data_prot.loc[i, "seq"]
                }
                row_rna = this.__fill_rna_row(rna_indices, dist_array, flag, err_message, transcript_id, this.__data_prot.loc[i, "seq"])
                results_dna.append(row_dna)
                results_rna.append(row_rna)
        df_dna = pd.DataFrame(results_dna)
        df_rna = pd.DataFrame(results_rna)
        this.__CreateDistanceFile(df_dna, splice_type,outputdir, basename)
        this.__CreateDistanceFile(df_rna, splice_type,outputdir, basename ,"RNA")

    def ComputeDistance(this, splice_type : str = "") -> dict[str : any]:
        """
        Method to gather all the distances computation in signe dictionnary
        """
        
        distance_calculator = {
            "A5SS_+": this.distanceA5SS,
            "A5SS_-": this.distanceA3SS,
            "A3SS_+": this.distanceA3SS,
            "A3SS_-": this.distanceA5SS,
            "RI_-": this.distanceRI ,
            "RI_+": this.distanceRI ,
            "SE_-": this.distanceSE ,
            "SE_+": this.distanceSE ,
            "MXE_-": this.distanceMXE,
            "MXE_+": this.distanceMXE,
            }

        return distance_calculator[splice_type](splice_type)
        
    def start(this, path_splicing : str, splice_type : str):
        """
        Method to start the computation of the distances
        path_prot : path to the file containing the protein data with the csv file at the end
        path_splicing : path to the folder containing the splicing data
        """
        this.__LoadDataSplicing(path_splicing)
        this.ComputeDistance(splice_type)




    def warmup_numba(this) -> None:
        # Jeu de données factice pour compiler Numba avant la parallélisation.
        coord = np.array([[100, 90]], dtype=np.int64)  # un seul couple
        exon_pos_list = [(50, 70)]
        _, _, _ = Distances.ComputeDistanceManual(coord, exon_pos_list)
        return

    
    def process_chunk(this,
                      df_chunk: pd.DataFrame,
                      df_splicing: pd.DataFrame,
                      comparison_couples: list[tuple[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Fonction qui applique la logique de start_manual sur un sous-ensemble (chunk) de df_ref.
        Elle renvoie deux DataFrame : (df_dna, df_rna)
        """
        results_dna = []
        results_rna = []

        # On parcourt le chunk local (au lieu de self.__data_prot complet)
        for i in range(len(df_chunk)):
            row_ref = df_chunk.iloc[i]
            
            # Récupération du Transcript et des exons
            transcript : pb.Transcript = this.bdd.transcript_by_id(row_ref["ensembl_id"])
            exon_pos_list = transcript.exon_intervals

            # Filtre sur df_splicing
            df_same_gene = df_splicing.loc[df_splicing["GeneID"] == row_ref["GeneID"]]

            for y in range(len(df_same_gene)):
                row_compare = df_same_gene.iloc[y]
                idx_couple = []
                for couple in comparison_couples:
                    array_coord = np.array([int(row_ref[couple[0]]), int(row_compare[couple[1]])])
                    idx_couple.append(array_coord)
                idx_couple = np.array(idx_couple)
                
                # Calcul des distances
                dist_array, flag_array, err_message_array = Distances.ComputeDistanceManual(idx_couple, exon_pos_list)

                # Construire la ligne "ADN"
                row_dna = {"transcript_ID": row_ref["ensembl_id"], "prot_seq": row_ref["seq"]}
                rna_indices = {}
                for i_couple, couple in enumerate(comparison_couples):
                    row_dna[f"{couple[0]}-{couple[1]}"] = dist_array[i_couple]
                    rna_indices[len(comparison_couples) + i_couple] = f"{couple[0]}-{couple[1]}"

                # Construire la ligne "ARN"
                row_rna = this.__fill_rna_row(
                    rna_indices,
                    dist_array,
                    flag_array,
                    err_message_array,
                    row_ref["ensembl_id"],
                    row_ref["seq"]
                )

                results_dna.append(row_dna)
                results_rna.append(row_rna)

        df_dna = pd.DataFrame(results_dna)
        df_rna = pd.DataFrame(results_rna)

        return df_dna, df_rna
    
    def parallel_start_manual(this,
                              df_ref: pd.DataFrame,
                              df_splicing: pd.DataFrame,
                              comparison_couples: list[tuple[str]],
                              output_dir : str,
                              output_basename : str = "manual",
                              n_cores : int = None
                             ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parallélise la logique de calcul sur df_ref en le découpant en plusieurs chunks.
        """
        if n_cores is None:
            n_cores = cpu_count()  # utilise tous les cœurs disponibles, ou définissez un nombre

        # Découper df_ref en n_cores (ou moins si df_ref est petit).
        df_ref = this.__FilterDataProt(df_ref)
        chunk_size = len(df_ref) // n_cores + 1
        chunks = [df_ref.iloc[i : i + chunk_size] for i in range(0, len(df_ref), chunk_size)]
        
        # Pour chaque chunk, on appelle process_chunk en parallèle
        with Pool(n_cores) as pool:
            results = pool.starmap(
                this.process_chunk,
                [(chunk, df_splicing, comparison_couples) for chunk in chunks]
            )

        # results est une liste de tuples (df_dna, df_rna) pour chaque chunk
        df_dna_concat = pd.concat([res[0] for res in results], ignore_index=True)
        df_rna_concat = pd.concat([res[1] for res in results], ignore_index=True)
        
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        df_rna_concat.to_csv(f"{output_dir}/rna_{output_basename}.csv", sep="\t", index=False)
        df_dna_concat.to_csv(f"{output_dir}/dna_{output_basename}.csv", sep="\t", index=False)

        return df_dna_concat, df_rna_concat


                 


if __name__ == "__main__":
    dist = Distances()
    # dist._LoadDataProt("data_filteredFMRP.tsv")
    # splice_types = [ "A5SS_+", "RI_+", "RI_-", "A5SS_-", "A3SS_+", "A3SS_-", "SE_-", "SE_+", "MXE_-", "MXE_+"]
    # for splice_type in splice_types:
    #     dist.start(path_splicing = "filteredRmats", splice_type = splice_type)
    #     print(f"Distances for {splice_type} computed")
    import time
    df_ref = pd.read_csv("data_filteredfinal.tsv", sep = "\t")
    df_2nd = pd.read_csv("filteredRmats/A5SS_+.csv", sep = "\t")
    couple_comparison = [("start_genomic", "shortSplice"), ("end_genomic", "longSplice"), ("start_genomic", "longSplice")]
    dist.warmup_numba() # compile numba avant la parallélisation
    debut = time.time()
    # dist.start_manual(df_ref, df_2nd, couple_comparison)
    dist.parallel_start_manual(df_ref, df_2nd, couple_comparison, n_cores = 4)
    print(time.time() - debut)
  