import ftplib
import os

def download_mouse_genome_ftp(output_folder):
    # Informations de connexion FTP
    ftp_server = "ftp.ncbi.nlm.nih.gov"
    ftp_path = "/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/"
    
    # Créer le dossier de sortie s'il n'existe pas
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Connexion au serveur FTP
    ftp = ftplib.FTP(ftp_server)
    ftp.login()
    ftp.cwd(ftp_path)
    
    # Lister les répertoires dans le chemin spécifié
    directories = ftp.nlst()
    
    for directory in directories:
        ftp.cwd(directory)
        files = ftp.nlst()
        
        for file in files:
            if file.endswith("_genomic.fna.gz"):  # Télécharger uniquement les fichiers .fna.gz
                local_filename = os.path.join(output_folder, file)
                with open(local_filename, 'wb') as local_file:
                    ftp.retrbinary(f"RETR {file}", local_file.write)
                print(f"Téléchargé : {file}")
        
        # Revenir au répertoire parent
        ftp.cwd("..")
    
    # Déconnexion du serveur FTP
    ftp.quit()
    print(f"Génome complet de la souris téléchargé dans le dossier {output_folder}")

# Télécharger le génome complet de la souris et le sauvegarder dans le dossier spécifié
download_mouse_genome_ftp("mouse_genome")