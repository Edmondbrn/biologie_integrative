# --------------------------------------------------------------------
# Utilise l'image Miniconda comme base (utile pour les dépendances python sinon faut installer anaconda depuis l'image ubuntu)
# --------------------------------------------------------------------
    FROM continuumio/miniconda3:latest

    ENV DEBIAN_FRONTEND=noninteractive
    
    # --------------------------------------------------------------------
    # Installe les dépendances bash nécessaires
    # --------------------------------------------------------------------
    RUN apt-get update && apt-get install -y \
        x11-apps \
        libgl1-mesa-glx \
        libegl1 \
        libdbus-1-3 \
        libx11-xcb1 \
        libxrender1 \
        libxi6 \
        libxtst6 \
        libxcomposite1 \
        libxcursor1 \
        libxdamage1 \
        libxfixes3 \
        libxrandr2 \
        libfontconfig1 \
        libfreetype6 \
        libxext6 \
        libx11-6 \
        libxkbcommon-x11-0 \
        libxcb1 \
        libxcb-util1 \
        libxcb-icccm4 \
        libxcb-image0 \
        libxcb-render-util0 \
        libxcb-randr0 \
        libxcb-shm0 \
        libxcb-xinerama0 \
        libxcb-keysyms1 \
        libxcb-xfixes0 \
        libxcb-shape0 \
        libxcb-sync1 \
        libxkbcommon0 \
        libsm6 \
        libxcb-cursor0 \
        && rm -rf /var/lib/apt/lists/*

    
    # --------------------------------------------------------------------
    # Création du répertoire app dans le conteneur et copie de l'environnement anconda
    # --------------------------------------------------------------------
    WORKDIR /app
    COPY env.yml /app
    
    # --------------------------------------------------------------------
    # Initialise l'environnement anaconda avec toutes les dépendances
    # --------------------------------------------------------------------
    RUN conda env create -f env.yml
    
    # --------------------------------------------------------------------
    # Active l'environnement anaconda
    # --------------------------------------------------------------------
    SHELL ["bash", "-c"]
    RUN echo "conda activate RepositionX" >> ~/.bashrc
    
    # --------------------------------------------------------------------
    # COpie tout le projet dans le conteneur
    # --------------------------------------------------------------------
    COPY . /app
    
    # --------------------------------------------------------------------
    # lance le projet
    # --------------------------------------------------------------------
    CMD ["bash", "-c", "source ~/.bashrc && python main.py"]
    