# MIC_Corrosion
This is the repository for the Dissertation "Nobel Microorganisms associated with MIC"
The first notebook is a visualisation the**Physicochemical** data analysis. These repository corresponds to the Session 3.3 Materials and Methods and 4.5 Results and Discussion in the main text.

| Notebook                    |Description                                                                                                  |
|----------------------------|--------------------------------------------------------------------------------------------------------------|
|**1_Splitdf.ipynb**     | Explaination of the data, then I separate the data into categories according to the failure study and the physicochemical study that distinguise three differnet degrees of corrosion.  The notebook analyse the type of data distribution. More details can be found in the main text of the dissertation (4.4.1). |
| **2_Filtering.ipynb**     | A notebook 1. Identify the most relevant bacteria genus(GIDs) amongts the 880 genus by filering out the percentage that is in a concentration of less than 2%, the genera are then sorted and an extracolumn tells the influence frecuency of those genera, lastly the data is divided into three traffic lights. A small statistic is done with the number of GID, count, Unique, top and frequency. |
| **3_PCA_RF_Feature.ipynb**      | This is a notebook with the analysis of the Principal Components of the variables of the data as well as the R F  . Refer to the main text (methods 3.3.5) and the results and discussion (Chapter 4.4.6) for details. |
| **4_Bacteria_Influencing_corrosion.ipynb**| Notebook for literature research on MIC, to differenciate the known and candidate bacteria on the selected_list of bacteria|
| **5_Sequences_find.ipynb**      | A notebook for invert search of the sequences from the curated data bacteria taxa which serve to make posterior a phylogenetic analysis and to know how this bacteria relate to each other from the evolutionary stand point.|
| **6_iTOL_tree.ipynb**      | A notebook for invert search of the sequences from the curated data bacteria taxa which serve to make posterior a phylogenetic analysis and to know how this bacteria relate to each other from the evolutionary stand point.|
| **7_Picrust_Functional.ipynb** | A notebook to analyse the functional and sequence relationships between newly identified bacteria and known corrosion-influencing microorganism, needs special environment and installation shown elsewhere*|

## Environment

Use the requirements file in this repo to create a new virtual environment for this task.

We have also added a [Makefile](Makefile) which has the recipe called 'setup' which will run all the commands for setting up the environment. Feel free to check and use if you are tired of copy pasting so many commands.

```BASH
make setup
#or
pyenv local 3.9.8
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
*Notebook __7_Picrust_Functional.ipynb__ needs to have installed miniconda, to be able to install
picrust2 here a short list of instructions: 
Install Miniconda:
Download Miniconda in home directory (~)
From Ubuntu terminal:
bashCopywget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
Accept license terms and default installation options

Create PICRUSt2 Environment:
bashCopyconda create -n picrust2 -c bioconda -c conda-forge picrust2
conda activate picrust2

Install Jupyter Support:
bashCopyconda install ipykernel
python -m ipykernel install --user --name=picrust2 --display-name="Python (picrust2)"

Install Required Packages:
bashCopyconda install -y pandas numpy biopython

Usage:
Navigate to repo: cd /path/to/repo
Activate environment: conda activate picrust2
Open VS Code: code .
In notebook: Select "Python (picrust2)" kernel

Note: PICRUSt2 environment is separate from the regular .venv environment.