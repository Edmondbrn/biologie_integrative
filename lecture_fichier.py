import pandas as pd
import requests
import json
from os import chdir, path

# Changer le répertoire de travail
chdir(path.dirname(__file__))

# Exemple de séquence et autres paramètres
sequence = "AGGATGATAGACGCGGTATTCACGAAAGACCCAGGCAGCAAGAAATGCACAAGCCTTTCCGAGGCTCAAATTT"
type_seq = "DNA"
database = "mm10"

# Construire l'URL
url = f"https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq={sequence}&type={type_seq}&db={database}&output=json"

# Envoyer la requête GET
response = requests.get(url)

# Vérifier si la requête a réussi
if response.status_code == 200:
    # Traiter la réponse JSON
    result = response.json()
    print(result)
    
    # Enregistrer la réponse JSON dans un fichier
    with open("result.json", "w") as json_file:
        json.dump(result, json_file, indent=4)
    
    # Récupérer le nom du gène du premier hit
    if "blat" in result and len(result["blat"]) > 0:
        first_hit = result["blat"][0]
        # Les champs sont définis dans "fields"
        fields = result["fields"]
        # Trouver l'index du champ "tName" (nom du gène)
        tname_index = fields.index("tName")
        gene_name = first_hit[tname_index]
        print(f"Nom du gène du premier hit : {gene_name}")
    else:
        print("Aucun hit trouvé dans la réponse.")
else:
    print(f"Erreur: {response.status_code}")