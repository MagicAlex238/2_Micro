# MIC_Corrosion
This is the repository for the Dissertation "Nobel Microorganisms associated with MIC"
The first notebook is a visualisation the**Physicochemical** data analysis. These repository corresponds to the Session 3.3 Materials and Methods and 4.5 Results and Discussion in the main text.


| Notebook                    |Description                                                                                                  |
|----------------------------|--------------------------------------------------------------------------------------------------------------|
|**1_Splitdf.ipynb**     | Explaination of the data, then I separate the data into categories according to the failure study and the physicochemical study that distinguise three differnet degrees of corrosion.  The notebook analyse the type of data distribution. More details can be found in the main text of the dissertation (4.4.1). |
| **2_Filtering.ipynb**     | A notebook 1. Identify the most relevant bacteria genus(GIDs) amongts the 880 genus by filering out the percentage that is in a concentration of less than 2%, the genera are then sorted and an extracolumn tells the influence frecuency of those genera, lastly the data is divided into three traffic lights. A small statistic is done with the number of GID, count, Unique, top and frequency. |
| **3_PCA_RF_Feature.ipynb**      | This is a notebook with the analysis of the Principal Components of the variables of the data as well as the R F  . Refer to the main text (methods 3.3.5) and the results and discussion (Chapter 4.4.6) for details. |
| **4_Bacteria_Influencing_corrosion.ipynb**| Notebook for literature research on MIC, to differenciate the known and candidate bacteria on the selected_list of bacteria, the notebook uses API calls to |
| **5_Sequences_qiime.ipynb**      | A notebook for invert search of the sequences from the curated data bacteria taxa which serve to make posterior a phylogenetic analysis and to know how this bacteria relate to each other from the evolutionary stand point|
| **6_Picrust_Functional.ipynb** | A notebook to analyse the functional and sequence relationships between newly identified bacteria and known corrosion-influencing microorganism, needs special environment and installation shown elsewhere*|


## Environment

my_mic_project/
├── .venv_general/           # For general notebooks (Notebooks 1, 2, 3 and 6)
├── .venv_notebook4/         # For Notebook 4 (Api calls with incompatible packages)
├── .venv_qiime/             # For Notebook 5 (Conda environment for QIIME)
├── notebooks/
│   ├── 1_splitdf.ipynb
│   ├── 2_Filtering.ipynb
│   ├── 3_Feature_selection.ipynb
│   ├── 4_Bacteria_Influencing_corrosion.ipynb
│   ├── 5_Sequences_qiime.ipynb
│   └── 6_picrust_functional.ipynb
├── requirements/
│   ├── general_requirements.txt
│   ├── notebook4_requirements.txt
│   ├── qiime_environment.yml  # For Conda
└── README.md
Following are instructions for local installation for vscode, done once. For Kaggle and Colab all installations have to be done every time the notebook is running.
## Instructions to install environment per notebook 1,2,3 and 6
1_ .venv_general/ 

```BASH

#For Linux
#
sudo apt update
sudo apt install python3.11 python3.11-venv python3.11-distutils
python3.11 -m venv .venv_general
source ./.venv_general/bin/activate
pip install --upgrade pip
pip install -r general_requirements.txt
python -m ipykernel install --user --name general_python11 --display-name "general_python11 (Python 3.11)"
```
3_ For powershell in Windows

```
python -m venv .venv_general
.\.venv_general\Scripts\Activate.ps1.
pip install --upgrade pip
pip install -r general_requirements.txt
python -m ipykernel install --user --name python11_general --display-name "python11_general"
```
## Instructions to install environment per notebook 4
1_ .venv_notebook4/ 

```BASH
2_For Linux
#
sudo apt update
sudo apt install python3.12 python3.12-venv python3.12-distutils
python3.12 -m venv .venv_notebook4
source .venv_notebook4/bin/activate
pip install --upgrade pip
pip install -r Notebook4_requirements.txt
python -m ipykernel install --user --name python12_N4 --display-name "python12_N4(Python 3.12)"

##  Instructions to install environments per .venv_qiime. Notebook 5_Sequences_qiime.ipynb`
qimme can only be used with linux installation, therefore the whole code has to be done under wsl2 linux for windows subsystem. It is necesary to install it first
# 1 Install WSL2 using Powershell
wsl --install # with administrations right
open new ubuntu and setup username and password
Install Miniconda within Ubuntu (WSL) terminal: Located on the new ubuntu terminal run
# 2 download the installation file for miniconda using bash shell
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -O miniconda_installer.sh
# 3 cloning the repository
git clone https://github.com/magicalex238/2_Micro.git # 
# 4 Navigate to the repository in wsl
cd 2_Micro
# 5. Download the requirements file for the qiime package using bash
wget https://data.qiime2.org/distro/core/qiime2-2023.7-py38-linux-conda.yml
# 6 Create the environment from the downloaded file in bash
conda env create -n qiime2-2023.7 -f qiime2-2023.7-py38-linux-conda.yml

# 7 Activate environment
conda activate qiime2-2023.7

# Verify installation
qiime --help

# 8 Install additional required packages
conda install -c bioconda -c conda-forge \
    biopython \
    pandas \
    numpy \
    matplotlib \
    seaborn

# 9 Import Visualisation tools
conda install -c bioconda -c conda-forge \
    itol-uploader \
    ete3

# 10 Install PICRUSt2 plugin
conda install -c conda-forge -c bioconda q2-picrust2

# 11 Install ipykernel  
conda install ipykernel  
# 12 Install the kernel for Jupyter 
python -m ipykernel install --user --name qiime2-2023.7 --display-name "Python (QIIME2)"

10 Usage:
Navigate to repo: cd /path/to/repo
Activate environment: conda activate qiime
Open VS Code: code .
In notebook: Select "Python (qiime)" kernel

