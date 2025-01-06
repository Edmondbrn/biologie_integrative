from numba import njit
import pandas as pd
import numpy as np
import pyensembl as pb
from multiprocessing import Pool, cpu_count
import os

ERROR_DICT = {4 : "Error while converting dna to rna",
                1 : "Not on the same transcript",
                2 : "Second coordinate not in the transcript",
                3 : "Error while searching for the second coordinates",
                0 : "No problemo"}

def FilterDataProt(df_prot : pd.DataFrame) -> pd.DataFrame:
    # enleve les lignes avec des str sur la colonne start_ensembl
    df_prot = df_prot.loc[df_prot["start_ensembl"].apply(lambda x: x.isnumeric())]
    return df_prot

def fill_rna_row(rna_indices : dict, dist_array : int, flag : bool, err_message :int , 
                    transcript_id : int, protein_sequence : str):
    """
    Méthode pour remplir le slignes des différents tableaux pour les distances ARN
    """
    row_rna = {}
    err_message = err_message[len(err_message)//2:]  # on ne garde que la 2e moitié qui correspond aux distances ARN
    for index, (i, col_name) in enumerate(rna_indices.items()):
        if err_message[index] != 0:
            # Remplacer la valeur par le message d’erreur ou un code
            row_rna[col_name] = f"ERROR_{ERROR_DICT[err_message[index]]}"
        elif flag[i]:
            # Ajouter un astérisque
            row_rna[col_name] = f"{dist_array[i]}*"
        else:
            # Valeur normale
            row_rna[col_name] = dist_array[i]
    row_rna["transcript_ID"]= transcript_id
    row_rna["prot_seq"]= protein_sequence
    return row_rna

@njit(cache = True, fastmath = True)
def get_intron_coord(exon_pos_list):
    """
    Retourne la liste des introns (start, end) à partir d'une liste
    d'exons [(ex1_start, ex1_end), (ex2_start, ex2_end), ...].
    Ex: si exon1 = [100, 120], exon2 = [140, 160], intron = [121, 139].
    """
    introns = []
    for i in range(len(exon_pos_list)-1):
        introns.append((exon_pos_list[i][1] + 1, exon_pos_list[i+1][0] - 1))
    return introns

@njit(cache = True, fastmath = True)
def get_intron_length(intron_start: int, intron_end: int) -> int:
    return max(0, intron_end - intron_start)  # évite négatif si chevauchement inattendu

@njit(cache = True, fastmath = True)
def _check_second_coordinate(max_coord: int,
                             i: int,
                             tot_exon: int,
                             dna_dist_abs: int,
                             rna_correction: int,
                             sign: int,
                             exon_pos_list: list[tuple[int, int]],
                             intron_pos_list: list[tuple[int, int]],
                             manual_flag: bool = False) -> tuple[int, bool, int]:
    """
    Sous-fonction qui gère la logique de 'max_coord' quand
    la coordonnée min est déjà positionnée (dans exon ou intron).
    max_cord : coordonnée de fin (ou de référence) du segment.
    i : index de l'exon où se trouve min_coord.
    tot_exon : nombre total d'exons.
    dna_dist_abs : distance absolue sur le génome.
    rna_correction : distance d'introns déjà retirée.
    sign : signe de la distance (1 ou -1).
    exon_pos_list : liste des exons.
    intron_pos_list : liste des introns.
    manual_flag : flag manuel (True) pour signaler un site dans un intron.
    
    Renvoie un tuple (distance_rna, has_star) ou une chaîne d'erreur.
    """
    # 1) Vérifier si max_coord est dans le même exon
    if exon_pos_list[i][0] <= max_coord-1 <= exon_pos_list[i][1]:
        # => Les 2 coordonnées sont dans le même exon
        bool_flag = manual_flag if manual_flag else False
        return sign * dna_dist_abs, bool_flag, 0

    # 2) Vérifier si c'est dans l'intron juste après l'exon i (si i existe)
    if i < len(intron_pos_list):
        if intron_pos_list[i][0] <= max_coord <= intron_pos_list[i][1]:
            # => La deuxième coordonnée est dans l'intron suivant
            return sign * dna_dist_abs, True , 0 # Signaler par un flag
    
    # 3) Sinon, on part du principe qu'on "saute" l'intron i pour aller au(x) suivant(s)
    # On retire la longueur de l'intron i (si on était dans l'exon i)
    if i < len(exon_pos_list) - 1:
        # Correction sur l'intron i
        rna_correction += get_intron_length(exon_pos_list[i][1], exon_pos_list[i+1][0])

    # 4) Boucle sur les exons suivants
    for y in range(i+1, tot_exon):
        # 4a) max_coord tombe dans l'exon y
        if exon_pos_list[y][0] <= max_coord <= exon_pos_list[y][1]:
            rna_dist_abs = dna_dist_abs - rna_correction
            bool_flag = manual_flag if manual_flag else False
            return sign * rna_dist_abs, bool_flag, 0
        
        # 4b) max_coord tombe dans l'intron suivant l'exon y
        if y < len(intron_pos_list):
            if intron_pos_list[y][0] <= max_coord <= intron_pos_list[y][1]:
                rna_dist_abs = dna_dist_abs - rna_correction
                return sign * rna_dist_abs, True, 0
        # 4c) Sinon, on ajoute la longueur de l'intron y pour aller à l'exon y+1
        if y < tot_exon - 1:
            rna_correction += get_intron_length(exon_pos_list[y][1], exon_pos_list[y+1][0])
        else:
            # Si on est déjà au dernier exon et on n'a pas trouvé "max_coord", c'est une erreur
            return 0, False, 2
    return 0, False, 3

@njit(cache = True, fastmath = True)
def convert_dna_to_rna( prot_coordinate: int,
                        splice_coordinate: int,
                        dna_distance: int,
                        exon_pos_list: list[tuple[int, int]]) -> tuple[int, bool, bool]:
    """
    Convertit la distance sur le génome (DNA) en distance sur l'ARN (RNA),
    en soustrayant les introns si nécessaire.
    Renvoie un tuple (distance_rna, has_star, code_erreur).
    
    has_star = True si on détecte un site dans un intron.
    """
    tot_exon = len(exon_pos_list)
    exon_pos_list = sorted(exon_pos_list)
    intron_pos_list = get_intron_coord(exon_pos_list)
    sign = 1 if dna_distance >= 0 else -1
    dna_dist_abs = abs(dna_distance)
    min_coord = min(prot_coordinate, splice_coordinate+1)
    max_coord = max(prot_coordinate, splice_coordinate+1)
    # Vérif : le splice_coordinate est-il dans la plage du transcript (ou exon/intron) ?
    # Ici, la condition "splice_coordinate < exon_pos_list[0][0]-1" ou ">" ... est la vôtre
    if splice_coordinate < (exon_pos_list[0][0] - 1) or splice_coordinate > (exon_pos_list[-1][1] + 1):
        return 0, False, 1
    # rna_correction accumule la longueur d'introns qu'on retire
    rna_correction = 0
    # On parcourt les exons pour voir où se situe min_coord
    for i in range(tot_exon):
        # CAS 1: min_coord est dans l'exon i
        if exon_pos_list[i][0] <= min_coord+1 <= exon_pos_list[i][1]:
            # Vérifie la deuxième coordonnée (max_coord)
            return _check_second_coordinate(max_coord, i, tot_exon,
                                            dna_dist_abs, rna_correction, sign,
                                            exon_pos_list, intron_pos_list)           
        # CAS 2: min_coord est dans l'intron i (si i existe dans intron_pos_list)
        if i < len(intron_pos_list):
            if intron_pos_list[i][0] <= min_coord <= intron_pos_list[i][1]:
                # 2a) Si max_coord est lui aussi dans le même intron i
                if intron_pos_list[i][0] <= max_coord <= intron_pos_list[i][1]:
                    # => Les deux coords dans le même intron => distance = dna_distance, star=True
                    return dna_distance, True, 0
                # Puis on check le 2e coord (max_coord) comme si on partait de exon i+1.
                return _check_second_coordinate(max_coord, i+1, tot_exon,
                                                dna_dist_abs, rna_correction, sign,
                                                exon_pos_list, intron_pos_list, True)
    # Si on sort de la boucle sans avoir rien retourné => problème
    return 0, False, 4


@njit(cache = True, fastmath = True)
def ComputeDistanceManual(coord : np.ndarray[np.ndarray[int]], 
                            exon_pos_list : list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # initialisation des tableaux numpy pour stocker les résultats
    nb_couple = coord.shape[0]
    dist = np.zeros(nb_couple * 2, dtype=np.int64)
    flag = np.zeros(nb_couple * 2, dtype=np.bool_)
    err_message = np.zeros(nb_couple * 2, dtype=np.int64)
    # parcours des différents couples de coordonnées pour les calculs
    for line in range(nb_couple):
        dist[line] = coord[line][0] - coord[line][1]
        dist[nb_couple + line], flag[nb_couple + line], err_message[nb_couple + line] = convert_dna_to_rna(
            coord[line][0],
            coord[line][1],
            dist[line],
            exon_pos_list
        )
    return dist, flag, err_message


def warmup_numba() -> None:
    # Jeu de données factice pour compiler Numba avant la parallélisation.
    coord = np.array([[100, 90]], dtype=np.int64)  # un seul couple
    exon_pos_list = [(50, 70)]
    _, _, _ = ComputeDistanceManual(coord, exon_pos_list)
    return


def process_chunk(df_chunk: pd.DataFrame,
                  df_splicing: pd.DataFrame,
                  comparison_couples: list[tuple[str]],
                  bdd : pb.EnsemblRelease
                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        transcript : pb.Transcript = bdd.transcript_by_id(row_ref["ensembl_id"])
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
            dist_array, flag_array, err_message_array = ComputeDistanceManual(idx_couple, exon_pos_list)
            # Construire la ligne "ADN"
            row_dna = {"transcript_ID": row_ref["ensembl_id"], "prot_seq": row_ref["seq"]}
            rna_indices = {}
            for i_couple, couple in enumerate(comparison_couples):
                row_dna[f"{couple[0]}-{couple[1]}"] = dist_array[i_couple]
                rna_indices[len(comparison_couples) + i_couple] = f"{couple[0]}-{couple[1]}"
            # Construire la ligne "ARN"
            row_rna = fill_rna_row(
                rna_indices,
                dist_array,
                flag_array,
                err_message_array,
                row_ref["ensembl_id"],
                row_ref["seq"])
            results_dna.append(row_dna)
            results_rna.append(row_rna)
    df_dna = pd.DataFrame(results_dna)
    df_rna = pd.DataFrame(results_rna)
    return df_dna, df_rna



def parallel_start_manual(df_ref: pd.DataFrame,
                          df_splicing: pd.DataFrame,
                          comparison_couples: list[tuple[str]],
                          bdd : pb.EnsemblRelease,
                          output_dir : str,
                          output_basename : str = "manual", 
                          n_cores : int = None,
                          progress_callback=None
                         ) -> None:
    """
    Parallélise la logique de calcul sur df_ref en le découpant en plusieurs chunks.
    Le paramètre progress_callback est une fonction de callback appelée
    à la fin de chaque chunk, pour mettre à jour la progression dans le GUI.
    """
    if n_cores is None:
        n_cores = cpu_count()
    warmup_numba() # Compilation pour les processus enfants
    df_ref = FilterDataProt(df_ref)
    chunk_size = len(df_ref) // n_cores + 1
    chunks = [df_ref.iloc[i : i + chunk_size] for i in range(0, len(df_ref), chunk_size)]
    # Pour accumuler les résultats finals
    results_dna = []
    results_rna = []
    def on_chunk_done(result):
        # Cette fonction est appelée dans le process principal 
        # dès qu'un chunk a fini de s’exécuter
        (df_dna_chunk, df_rna_chunk) = result
        results_dna.append(df_dna_chunk)
        results_rna.append(df_rna_chunk)
        # On peut appeler la callback pour la progression
        if progress_callback is not None:
            progress_callback(len(df_dna_chunk)) 
            # on fait progresser de +NbLignesChunk 
    with Pool(n_cores) as pool:
        # Lancement asynchrone de chaque chunk
        for chunk in chunks:
            pool.apply_async(
                process_chunk,
                args=(chunk, df_splicing, comparison_couples, bdd),
                callback=on_chunk_done
            )
        # On bloque jusqu'à ce que tout soit fini
        pool.close()
        pool.join()
    # results_dna / results_rna sont maintenant alimentés par le callback
    df_dna_concat = pd.concat(results_dna, ignore_index=True)
    df_rna_concat = pd.concat(results_rna, ignore_index=True)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    df_rna_concat.to_csv(f"{output_dir}/rna_{output_basename}.csv", sep="\t", index=False)
    df_dna_concat.to_csv(f"{output_dir}/dna_{output_basename}.csv", sep="\t", index=False)
    return



def process_chunk_splicing(df_prot: pd.DataFrame,
                            df_splicing: pd.DataFrame,
                            comparison_couples: list[tuple[str, str]],
                            bdd : pb.EnsemblRelease,
                            splice_type : str
                            ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fonction qui applique la logique de start_manual sur un sous-ensemble (chunk) de df_ref.
    Elle renvoie deux DataFrame : (df_dna, df_rna)
    """
    results_dna = []
    results_rna = []
    # On parcourt data prot complet ici
    for i in range(len(df_prot)):
        row_ref = df_prot.iloc[i]
        # Récupération du Transcript et des exons
        transcript : pb.Transcript = bdd.transcript_by_id(row_ref["ensembl_id"])
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
            dist_array, flag_array, err_message_array = ComputeDistanceManual(idx_couple, exon_pos_list)
            # Construire la ligne "ADN"
            row_dna = {"transcript_ID": row_ref["ensembl_id"], "prot_seq": row_ref["seq"]}
            rna_indices = {}
            for i_couple, couple in enumerate(comparison_couples):
                row_dna[f"{couple[0]}-{couple[1]}"] = dist_array[i_couple]
                rna_indices[len(comparison_couples) + i_couple] = f"{couple[0]}-{couple[1]}"
            # Construire la ligne "ARN"
            row_rna = fill_rna_row(
                rna_indices,
                dist_array,
                flag_array,
                err_message_array,
                row_ref["ensembl_id"],
                row_ref["seq"])
            results_dna.append(row_dna)
            results_rna.append(row_rna)
    df_dna = pd.DataFrame(results_dna)
    df_rna = pd.DataFrame(results_rna)
    return df_dna, df_rna, splice_type



def parallel_start_manual_all(df_ref: pd.DataFrame,
                            input_dfs: dict[str : pd.DataFrame],
                            comparison_couples: dict[str : list[tuple[str, str]]],
                            bdd : pb.EnsemblRelease,
                            output_dir : str,
                            output_basename : str = "manual", 
                            n_cores : int = None,
                            progress_callback=None
                            ) -> None:
    """
    Parallélise la logique de calcul sur df_ref en le découpant en plusieurs chunks.
    Le paramètre progress_callback est une fonction de callback appelée
    à la fin de chaque chunk, pour mettre à jour la progression dans le GUI.
    """
    if n_cores is None:
        n_cores = len(input_dfs)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    warmup_numba() # Compilation pour les processus enfants
    df_ref = FilterDataProt(df_ref)
    # chunk_size = len(df_ref) // n_cores + 1
    # chunks = [df_ref.iloc[i : i + chunk_size] for i in range(0, len(df_ref), chunk_size)]

    
    def on_chunk_done(result):
        # Cette fonction est appelée dans le process principal 
        # dès qu'un chunk a fini de s’exécuter
        (df_dna_chunk, df_rna_chunk, splice_type) = result
        df_dna_chunk.to_csv(f"{output_dir}/dna_{splice_type}_{output_basename}.csv", sep="\t", index=False)
        df_rna_chunk.to_csv(f"{output_dir}/rna_{splice_type}_{output_basename}.csv", sep="\t", index=False)

        # On peut appeler la callback pour la progression
        if progress_callback is not None:
            progress_callback(len(df_dna_chunk)) 
            # on fait progresser de +NbLignesChunk 
    with Pool(n_cores) as pool:
        # Lancement asynchrone de chaque chunk
        for splice_type, splicing_chunk in input_dfs.items():
            pool.apply_async(
                process_chunk_splicing,
                args=(df_ref, splicing_chunk, comparison_couples[splice_type], bdd, splice_type),
                callback=on_chunk_done
            )
        # On bloque jusqu'à ce que tout soit fini
        pool.close()
        pool.join()
    return




if __name__ == "__main__":
    # bdd = pb.EnsemblRelease(species="mus_musculus", release=102)
    exons = [(100, 120), (140,160)]
    data_test = [
        # -- A. Les deux coordonnées dans le même exon --
        # 1) Les deux dans Exon1 [100..120]
        (105, 110, "Les deux coordonnées dans l'exon1 (segment interne)"),
        (100, 120, "Les deux coordonnées couvrent entièrement exon1"),

        # 2) Les deux dans Exon2 [140..160]
        (145, 150, "Les deux coordonnées dans l'exon2 (segment interne)"),
        (140, 160, "Les deux coordonnées couvrent entièrement exon2"),

        # -- B. L’une dans Exon1, l’autre dans Exon2 => traverse intron complet --
        (120, 140, "Juste à la frontière Exon1/Intron1 et Exon2"),
        (100, 160, "De l'extrémité de l'exon1 à l'extrémité exon2 (traverse intron complet)"),

        # -- C. Les deux dans l’intron [121..139] --
        (125, 130, "Les deux coordonnées sont dans l'intron (segment court)"),
        (121, 139, "Segment couvre entièrement intron1"),

        # -- D. Une coordonnée dans l’intron, l’autre dans un exon --
        (115, 130, "ref dans exon1, end dans intron1"),
        (130, 150, "ref dans intron1, end dans exon2"),
        (120, 130, "ref tout juste fin exon1, end dans intron1"),
        (139, 140, "ref tout juste fin intron1, end tout début exon2"),

        # -- E. Cas "extrêmes" / en dehors du gène --
        (90, 95,  "Les deux coordonnées avant exon1 => pas de traversée intron"),
        (165, 170, "Les deux coordonnées après exon2 => pas de traversée intron"),
        (95, 145, "Traverse partiellement exon1 et l'intron pour arriver dans exon2"),
    ]
    for i in range(0, len(data_test)):
        ref_coord = data_test[i][0]
        end_corrd = data_test[i][1]
        dna_distance = ref_coord - end_corrd
        result = convert_dna_to_rna(ref_coord, end_corrd, dna_distance, exons)
        print(f"TEST : {data_test[i][2]}\n\
                ref_coordinate = {ref_coord} end_coordinate : {end_corrd}\n\
                dna_distance = {dna_distance}\n\
                Résultat : {result}\n\
                Exons : {exons}")

  