from numba import njit

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
    elif end_coordinate >= exon_pos[0] and end_coordinate <= exon_pos[1]:
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
def convert_dna_to_rna(ref_coordinate: int, end_coordinate: int, dna_distance: int, exon_pos_list : list[tuple[int, int]]) -> int:
    """
    Convertit la distance sur le génome (DNA) en distance sur l'ARN (RNA),
    en soustrayant les introns si nécessaire.
    """

    for i in range(len(exon_pos_list)):
        exon = exon_pos_list[i]
        if is_both_between(ref_coordinate, end_coordinate, exon):
            return dna_distance
        else:
            check, unmatch_coord = is_only_one_between(ref_coordinate, end_coordinate, exon)
            if check:
                # On détermine l'intron suivant
                correction = 0
                for j in range(i+1, len(exon_pos_list)):
                    exon_j = exon_pos_list[j]
                    if is_in_exon(unmatch_coord, exon_j):
                        correction += get_intron_length(exon_pos_list[j-1][1], exon_j[0])
                        dna_distance = add_or_subtract_intron_length(dna_distance, correction)
                        return dna_distance
                    elif is_in_previous_intron(unmatch_coord, exon_j):
                        correction += get_intron_length(exon_pos_list[j-1][1], unmatch_coord)
                        dna_distance = add_or_subtract_intron_length(dna_distance, correction)
                        return dna_distance
                    elif j == len(exon_pos_list) - 1:
                        dna_distance = add_or_subtract_intron_length(dna_distance, correction)
                        return dna_distance
                    else:
                        correction += get_intron_length(exon_pos_list[j-1][1], exon_pos_list[j][0])
    return dna_distance  # Par défaut si aucune condition ne s'applique

if __name__ == "__main__":
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
