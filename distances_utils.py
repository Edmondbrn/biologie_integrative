from numba import njit
import pyensembl  as pb
from numpy import Inf
import numpy as np

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
    intron_list.append((exon_pos_list[-1][1]+1, Inf))
    return intron_list

@njit(fastmath=True)
def convert_dna_to_rna(ref_coordinate: int, end_coordinate: int, dna_distance: int, exon_pos_list : list[tuple[int, int]]) -> int:
    """
    Convertit la distance sur le génome (DNA) en distance sur l'ARN (RNA),
    en soustrayant les introns si nécessaire.
    """
    exon_pos_list = sorted(exon_pos_list)
    intron_pos_list = get_intron_coord(exon_pos_list)
    min_coord = min(ref_coordinate, end_coordinate)
    max_coord = max(ref_coordinate, end_coordinate)
    rna_correction = 0
    for i in range(len(exon_pos_list)):
        exon = exon_pos_list[i]
        intron = intron_pos_list[i]
        if is_both_between(ref_coordinate, end_coordinate, exon):
            return dna_distance
        elif is_in_exon(min_coord, exon) or is_in_intron(min_coord, intron):
            if is_in_intron(max_coord, intron):
                rna_correction += get_intron_length(exon[1], max_coord)
                return add_or_subtract_intron_length(dna_distance, rna_correction)
            elif is_in_exon(max_coord, exon):
                return dna_distance
            else:
                rna_correction = min(intron[1], max_coord) - max(min_coord, exon[1])
                for j in range(i+1, len(exon_pos_list)):
                    exon_pos_j = exon_pos_list[j]
                    intron_pos_j = intron_pos_list[j]
                    if is_in_exon(max_coord, exon_pos_j):
                        return add_or_subtract_intron_length(dna_distance, rna_correction)
                    elif is_in_intron(max_coord, intron_pos_j):
                        return add_or_subtract_intron_length(dna_distance, rna_correction)
                    else:
                        rna_correction += get_intron_length(exon_pos_j[1], intron_pos_j[1]+1)
    return 1000000000
                        
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

if __name__ == "__main__":
    bdd = pb.EnsemblRelease(species="mus_musculus", release=102)
    print("Test 1")
    exon_pos_list = [(0, 10), (15, 20), (25, 30)]
    ref_coordinate = 5
    end_coordinate = 18
    dna_distance = -13
    print(f"Exepected RNA distance: -8")
    print(convert_dna_to_rna(ref_coordinate, end_coordinate, dna_distance, exon_pos_list))
    print("="*80)
    print("Test 2") # si les deux positions sont dans le même exon
    exon_pos_list = [(0, 10), (15, 20), (25, 30)]
    ref_coordinate = 5
    end_coordinate = 8
    dna_distance = -3
    print(f"Exepected RNA distance: -3")
    print(convert_dna_to_rna(ref_coordinate, end_coordinate, dna_distance, exon_pos_list))
    print("="*80)
    print("Test 3") # si la position de fin est avant la position de début
    exon_pos_list = [(0, 10), (15, 20), (25, 30)]
    ref_coordinate = 16
    end_coordinate = 3
    dna_distance = 13
    print(f"Exepected RNA distance: 13")
    print(convert_dna_to_rna(ref_coordinate, end_coordinate, dna_distance, exon_pos_list))
    print("="*80)
    print("Test 4") # s'il ya plus d'un exon d'écart
    exon_pos_list = [(0, 10), (15, 20), (25, 30)]
    ref_coordinate = 5
    end_coordinate = 28
    dna_distance = -23
    print(f"Exepected RNA distance: -13")
    print(convert_dna_to_rna(ref_coordinate, end_coordinate, dna_distance, exon_pos_list))
    print("Debug RI coordinate")
    transcript_id = "ENSMUST00000110168"
    coord1 = 58069393
    coord2 = 58052223
    dna_distance = coord1 - coord2
    transcript : pb.Transcript = bdd.transcript_by_id(transcript_id)
    exon_pos = transcript.exon_intervals
    print(f"1st coordinate: {coord1}")
    print(f"2nd coordinate: {coord2}")
    print(f"Transcript ID: {transcript_id}")
    print(f"DNA distance: {dna_distance}")
    print(f"RNA distance: {convert_dna_to_rna(coord1, coord2, dna_distance, exon_pos)}")
    # print("test new function")
    # convert_dna_to_rna(coord1, coord2, dna_distance, exon_pos)
