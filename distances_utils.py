from numba import njit


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
    if exon_pos_list[i][0] <= max_coord <= exon_pos_list[i][1]:
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
        if exon_pos_list[i][0] <= min_coord <= exon_pos_list[i][1]:
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
                
                # Ajout correction intron i :
                # if (i+1) < tot_exon:
                    # rna_correction += get_intron_length(min_coord, intron_pos_list[i][1])
                
                # Puis on check le 2e coord (max_coord) comme si on partait de exon i+1.
                return _check_second_coordinate(max_coord, i+1, tot_exon,
                                                dna_dist_abs, rna_correction, sign,
                                                exon_pos_list, intron_pos_list, True)
    # Si on sort de la boucle sans avoir rien retourné => problème
    return 0, False, 4





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

  