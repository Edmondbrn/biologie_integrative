
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def create_file(id_list : list[str], file_path : str):
    liste_size = len(id_list)
    if liste_size == 0:
        return
    elif liste_size > 500:
        print("The list is too long. It will be cut to 500 elements.")
        for i in range(liste_size // 500 +1):
            print(liste_size // 500)
            split_id = id_list[i*500:(i+1)*500]
            write_file(split_id, file_path + "_" + str(i))
    else:
        write_file(id_list, file_path)
    
def write_file(id_list : list[str], file_path : str):
    with open(file_path, 'w') as file:
        for id in id_list:
            id = id.split(".")[0]
            file.write(id + "\n")

def create_id_dict(file_path : str) -> dict[str, str]:
    id_dict = {"NM": [], "NR": [], "XM": [], "XR": []}
    with open(file_path, 'r') as file:
        for line in file:
            id = line.strip()
            id_dict[id.split("_")[0]].append(id)
    return id_dict

path = "id_ncbi.txt"
dict_id = create_id_dict(path)
for key in dict_id:
    create_file(dict_id[key], "output/" + key)