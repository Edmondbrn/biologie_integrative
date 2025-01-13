
# RepostionX
<img src="src/Ressources/Icones/fugue-icons-3.5.6/icons/BI_logo.png" alt="BI Logo" width="200">

RepositionX is a `python3` free and opensource software wich allow you to compute distances between two sites on RNA or DNA. RepositionX is a final student year project at Polytech Nice-Sophia bio-engineering school in partnership with the labaratory IPMC memeber of the CNRS.

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
[![python](https://img.shields.io/badge/Python-3.10.15-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)

![Anaconda](https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white)
![Static Badge](https://img.shields.io/badge/Anaconda-v24.5.0-brightgreen)

![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![Static Badge](https://img.shields.io/badge/Docker-v27.4.1-brightgreen?logoColor=(0%2C0%2C255))
## Installation

```batch
    python --version
```
Or
```bash
    python3 --version
```
    

- Windows: 
There is a .exe file which launch the project. But you can also use the .bat file which will install all the python dependencies and run the program.


- Linux:
The project contains a Anaconda environment file. You can install the virtual environement with th command :

```bash
  conda env create -f env.yml
```

And launch the project with :

```bash
  conda activate <env_name>
  python3 main.py
```

- Docker:
A docker file is also furnished with the project. Linux users can directly build the image with :

```bash
docker build -t <your-username>/RepositionX dockerfile
docker image ls
```

After the building you can run the docker with:

```bash
bash boot_docker.sh
```

Note that if you want to get back the ouput file of the program, you wille have to mount your output repository inside the docker. See the commented line inside boot_docker.


    
## Features

- Cross platform
- Convert NCBI transcript IDs into Ensembl IDs
- Convert mRNA coordinates into DNA coordinates
- Filter and restructure Rmats ouput files
- Compute DNA distances and mRNA distances
- Visualise protein fixation sites ont raw RNA


## Documentation

This part will be dedicated to explain how to use the software.

- Home view
<img src="src/Ressources/readme_data/home.png" alt="home" width="400">

When you start the program, a window like this will open. From this window you can load datasets from the menu `File`, or from the table icon. If you don't want to load a datasets directly, you can click on the `Actions` menu to get access to the features.


- Compute Distances

This menu will trigger various sub-menus to let you choose what type of calculation you want to. You can choose between :

    1. A5SS alternative splicing sites
    2. A3SS alternative splicing sites
    3. RI alternative splicing sites
    4. MXE alternative splicing sites
    5. SE alternative splicing sites
    6. All types of alternative splicing sites
    7. Manual mode

After your choice, a window of this style will open:

<img src="src/Ressources/readme_data/distance_menu1.png" alt="home" width="400">

You have a section to load the reference datasets. The calculations will be based on these coordinates (ref_coordinates - 2nd_coordinates). Also for the comparaison datasets. You can load excel and CSV files. Note that for CSV files, we will scan the file to set the separator characater. If it fails, please select yourself the character.

**Datasets have to respect some rules:**
- Reference dataset has to contain the `GeneID`column with Ensembl gene identifyers for the element and `ensembl_id` column which contains the Ensembl transcript identifyer.
- Comparison dataset has to contain the `GeneID`column with Ensembl gene identifyers for the element


After that you have to select an output directory ans an output file name. If nothing is given, the output files will be in the root prooject with just `dna`or `rna`prefix. If everythin is done, you can click and `Validate`.

<img src="src/Ressources/readme_data/distance_menu2.png" alt="home" width="600">

Some elements will be added to the window. Now you can select the columns to compare with the combo boxes. If you have chosen a specific splicing type or all, the comparison coluns will be automaticaly added.

**Please note that you have to specify a column with __only__ numeric values. Remove strings or ambiguous values. Otherwise, you will have an error message.**

You can add all the combinaison you want, but not twice the same (Warning message if this case). If you have done a mistake, you can remove the last added combinaison with the red button. Finally you can launch the computation with the `Compare` button. A progress bar will be displayed to warn you about the calculation progress. A little moment can occured before the starting of calculation because of compilation of some functions.

**Note that the calculations can be parallelised by activating the multiprocessing. A spinbox will appeare to ask you how many processes you want. This is highly not recommended an the .exe version.**
## Authors

**First author**
- [@Berne Edmond](https://github.com/Edmondbrn)

**Other authors**
- [@Martin François](https://www.github.com/octokatherine)
- [@Rouget Simon](https://www.github.com/octokatherine)


