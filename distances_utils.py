from numba import njit
import pyensembl  as pb
from numpy import Inf
import pytest

@njit(fastmath=True)
def is_both_between(start_coordinate: int, end_coordinate: int, exon_pos: tuple) -> bool:
    """
    Vérifie si deux positions sont dans le même exon.
    """
    return (
        start_coordinate >= exon_pos[0] and start_coordinate <= exon_pos[1]
    ) and (
        end_coordinate >= exon_pos[0] and end_coordinate <= exon_pos[1]
    )

@njit(fastmath=True)
def is_only_one_between(start_coordinate: int, end_coordinate: int, exon_pos: tuple):
    """
    Vérifie si une seule des deux positions se trouve dans un exon.
    Retourne (True, autre_coord) si c'est le cas, sinon (False, None).
    """
    if start_coordinate >= exon_pos[0] and start_coordinate <= exon_pos[1]:
        return (True, end_coordinate)
    elif end_coordinate >= exon_pos[0] and end_coordinate -1 <= exon_pos[1]:
        return (True, start_coordinate)
    else:
        return (False, None)
    
@njit(fastmath=True)
def is_in_exon(coordinate: int, exon_pos: tuple) -> bool:
    """
    Vérifie si une position se trouve dans un exon.
    """
    return exon_pos[0] <= coordinate <= exon_pos[1]

@njit(fastmath=True)
def is_in_intron(coordinate: int, intron_pos: tuple) -> bool:
    """
    Vérifie si une position se trouve dans un exon.
    """
    return intron_pos[0] <= coordinate <= intron_pos[1]

@njit(fastmath=True)
def get_intron_length(exon1_end: int, exon2_start: int) -> int:
    """
    Retourne la longueur de l'intron entre deux exons.
    """
    return exon2_start - exon1_end

@njit(fastmath=True)
def add_or_subtract_intron_length(dna_distance: int, intron_length : int) -> int:
    """
    Ajoute ou soustrait la longueur des introns à la distance sur le génome (DNA).
    """
    if dna_distance < 0:
        return dna_distance + intron_length
    else:
        return dna_distance - intron_length

@njit(fastmath=True)
def is_in_previous_intron(coordinate: int, exon_pos: tuple) -> bool:
    """
    Vérifie si une position se trouve dans l'intron précédent un exon.
    """
    return coordinate < exon_pos[0]

@njit(fastmath=True)
def get_intron_coord(exon_pos_list : tuple[list[int, int]])-> tuple[int, int]:
    """
    Retourne les coordonnées de l'intron entre deux exons.
    """
    intron_list = list()
    for j in range(len(exon_pos_list)-1):
        intron_list.append((exon_pos_list[j][1]+1, exon_pos_list[j+1][0]-1))
    # intron_list.append((exon_pos_list[-1][1]+1, Inf))
    return intron_list

# @njit(fastmath=True)
# def convert_dna_to_rna(ref_coordinate: int, end_coordinate: int, dna_distance: int, exon_pos_list : list[tuple[int, int]]) -> int:
#     """
#     Convertit la distance sur le génome (DNA) en distance sur l'ARN (RNA),
#     en soustrayant les introns si nécessaire.
#     """
#     exon_pos_list = sorted(exon_pos_list)
#     intron_pos_list = get_intron_coord(exon_pos_list)
#     min_coord = min(ref_coordinate, end_coordinate)
#     max_coord = max(ref_coordinate, end_coordinate)
#     rna_correction = 0
#     for i in range(len(exon_pos_list)):
#         exon = exon_pos_list[i]
#         intron = intron_pos_list[i]
#         if is_both_between(ref_coordinate, end_coordinate, exon):
#             return dna_distance
#         elif is_in_exon(min_coord, exon) or is_in_intron(min_coord, intron):
#             if is_in_intron(max_coord, intron):
#                 rna_correction += get_intron_length(exon[1], max_coord)
#                 return add_or_subtract_intron_length(dna_distance, rna_correction)
#             elif is_in_exon(max_coord, exon):
#                 return dna_distance
#             else:
#                 rna_correction = min(intron[1], max_coord) - max(min_coord, exon[1])
#                 for j in range(i+1, len(exon_pos_list)):
#                     exon_pos_j = exon_pos_list[j]
#                     intron_pos_j = intron_pos_list[j]
#                     if is_in_exon(max_coord, exon_pos_j):
#                         return add_or_subtract_intron_length(dna_distance, rna_correction)
#                     elif is_in_intron(max_coord, intron_pos_j):
#                         return add_or_subtract_intron_length(dna_distance, rna_correction)
#                     else:
#                         rna_correction += get_intron_length(exon_pos_j[1], intron_pos_j[1]+1)
#     return 1000000000
                        
# @njit(fastmath=True)
# def convert_dna_to_rna(ref_coordinate: int, end_coordinate: int, dna_distance: int, exon_pos_list : list[tuple[int, int]]) -> int:
#     """
#     Convertit la distance sur le génome (DNA) en distance sur l'ARN (RNA),
#     en soustrayant les introns si nécessaire.
#     """
#     exon_pos_list = sorted(exon_pos_list)
#     for i in range(len(exon_pos_list)):
#         exon = exon_pos_list[i]
#         if is_both_between(ref_coordinate, end_coordinate +1, exon): # si les 2 coordonnées sont sur le même exon
#             return dna_distance
#         else:
#             # print(f"Ref coord: {ref_coordinate}, End coord: {end_coordinate}, Exon: {exon[0]} : {exon[1]}")
#             check, unmatch_coord = is_only_one_between(ref_coordinate, end_coordinate+1, exon) # si une des deux est sur l'exon étudié
#             if check:
#                 # On détermine l'intron suivant
#                 correction = 0
#                 # print(f"Exon pos start: {exon[0]} : {exon[1]}")
#                 for j in range(i+1, len(exon_pos_list)): # parcours des exons suivants
#                     # print(f"Exon pos : {exon_pos_list[j][0]} : {exon_pos_list[j][1]}")
#                     # print(f"Correction: {correction}")
#                     exon_j = exon_pos_list[j]
#                     if is_in_exon(unmatch_coord, exon_j): # si l'autre coordonnée est dans l'exon suivant
#                         # print("Test is in exon")
#                         correction += get_intron_length(exon_pos_list[j-1][1], exon_j[0])
#                         dna_distance = add_or_subtract_intron_length(dna_distance, correction)
#                         return dna_distance
#                     elif is_in_previous_intron(unmatch_coord, exon_j): # si l'autre coordonnée est dans l'intron précédent de l'exon suivant
#                         # print("Test is in previous intron")
#                         correction += get_intron_length(exon_pos_list[j-1][1], unmatch_coord)
#                         # print(f"Correction 2: {correction}")
#                         dna_distance = add_or_subtract_intron_length(dna_distance, correction)
#                         return dna_distance
#                     elif j == len(exon_pos_list) - 1: # si on est à la fin de la liste des exons
#                         # print("Test is at the end")
#                         dna_distance = add_or_subtract_intron_length(dna_distance, correction)
#                         return dna_distance
#                     else: #~ si l'autre coordonnée n'est pas dans l'exon suivant
#                         # print("Test is not in exon")
#                         correction += get_intron_length(exon_pos_list[j-1][1], exon_pos_list[j][0])
#     return dna_distance  # Par défaut si aucune condition ne s'applique

@njit(fastmath=True)
def verify_count(dna_distance, rna_distance):
    print(f"Error, rna : {rna_distance}, dna : {dna_distance}") if abs(rna_distance) > abs(dna_distance) else None




# def convert_dna_to_rna(dna_pos : int, transcript : pb.Transcript) -> int :
#     rna_pos = 0
#     exons = transcript.exon_intervals

#     for exon_start, exon_end in exons:
#         print(exon_start)
#         if exon_start <= dna_pos <= exon_end:
#             offset = dna_pos - exon_start
#             return rna_pos + offset
#         elif dna_pos > exon_end:
#             rna_pos += exon_end - exon_start
#         else:
#             return "Error dna position is before than first exon start site"


# @njit(fastmath=True)
# def convert_dna_to_rna(prot_coordinate: int, splice_coordinate: int, dna_distance: int, exon_pos_list : list[tuple[int, int]]) -> int:
#     """
#     Convertit la distance sur le génome (DNA) en distance sur l'ARN (RNA),
#     en soustrayant les introns si nécessaire.
#     """
#     tot_exon = len(exon_pos_list)
#     exon_pos_list = sorted(exon_pos_list)
#     intron_pos_list = get_intron_coord(exon_pos_list)
#     rna_correction = 0
#     sign = 1 if dna_distance >= 0 else -1
#     dna_dist_abs = abs(dna_distance)
#     min_coord = min(prot_coordinate, splice_coordinate)
#     max_coord = max(prot_coordinate, splice_coordinate)
#     if splice_coordinate < exon_pos_list[0][0]-1 or splice_coordinate > exon_pos_list[-1][1]+1:
#         return "Error splice site is not on the same transcript as protein fixation site"
#     else:
#         for i in range(tot_exon):
#             if exon_pos_list[i][0] <= min_coord <= exon_pos_list[i][1]: # si la coordonnées la plus petite est dans l'exon
#                 if exon_pos_list[i][0] <= max_coord <= exon_pos_list[i][1]: # si la plus grande l'est aussi
#                     return dna_distance, False # on ne fait rien à la distance ADN alors
#                 elif  i != tot_exon and intron_pos_list[i][0] <= max_coord <= intron_pos_list[i][1] : # si le deuxième site est dans l'intron suivant
#                     return str(dna_distance) , True # si le deuxième site est dans un intron, on garder la distance ADN mais on met une * pour avertir l'utilisateur
#                 else : # si le deuxième site est dans une structure plus loin
#                     rna_correction += get_intron_length(exon_pos_list[i][1], exon_pos_list[i+1][0])
#                     for y in range(i+1, tot_exon):
#                         if exon_pos_list[y][0] <= max_coord <= exon_pos_list[y][1]: # si la deuxième coordonnées est dans l'exon suivant
#                             rna_dist_abs = dna_dist_abs - rna_correction # on retire  la distance des introns passés
#                             return sign * rna_dist_abs, False
#                         elif y != tot_exon and intron_pos_list[y][0] <= max_coord <= intron_pos_list[y][1]:
#                             rna_dist_abs = dna_dist_abs - rna_correction
#                             return sign * rna_dist_abs, True # allerte l'utilisateur si la deuxième coordonnées est dans un intron
#                         else:
#                             rna_correction += get_intron_length(exon_pos_list[y][1], exon_pos_list[y+1][0])
#                     break
#             elif intron_pos_list[i][0] <= min_coord <= intron_pos_list[i][1]:
#                 pass
#             else:
#                 continue
#         return "Error while computing RNA distance"


@njit
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

@njit
def get_intron_length(intron_start: int, intron_end: int) -> int:
    return max(0, intron_end - intron_start)  # évite négatif si chevauchement inattendu

@njit
def _check_second_coordinate(
    max_coord: int,
    i: int,
    tot_exon: int,
    dna_dist_abs: int,
    rna_correction: int,
    sign: int,
    exon_pos_list: list[tuple[int, int]],
    intron_pos_list: list[tuple[int, int]],
    manual_flag: bool = False
):
    """
    Sous-fonction qui gère la logique de 'max_coord' quand
    la coordonnée min est déjà positionnée (dans exon ou intron).
    
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

@njit
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
                if (i+1) < tot_exon:
                    rna_correction += get_intron_length(min_coord, intron_pos_list[i][1])
                
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

  