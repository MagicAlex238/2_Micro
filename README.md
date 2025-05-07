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
| **7_visual_protein.ipynb** | A notebook to visualise the computational results of the picrust2 prediction|

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
*Notebook Instructions for

# 5_Sequences_qiime.ipynb Notebook
# download the installation file
wget https://data.qiime2.org/distro/core/qiime2-2023.7-py38-linux-conda.yml

# Create the environment from the downloaded file
conda env create -n qiime2-2023.7 -f qiime2-2023.7-py38-linux-conda.yml

# Create and activate environment
conda create -n qiime2-2023.7 -c conda-forge -c bioconda qiime2
conda activate qiime2-2023.7

# Verify installation
qiime --help

# Install additional required packages
conda install -c bioconda -c conda-forge \
    biopython \
    pandas \
    numpy \
    matplotlib \
    seaborn

# Import Visualisation tools
conda install -c bioconda -c conda-forge \
    itol-uploader \
    ete3

# Install PICRUSt2 plugin
conda install -c conda-forge -c bioconda q2-picrust2

# Install ipykernel  
conda install ipykernel  
# Install the kernel for Jupyter 
python -m ipykernel install --user --name qiime2-2023.7 --display-name "Python (QIIME2)"

Usage:
Navigate to repo: cd /path/to/repo
Activate environment: conda activate qiime
Open VS Code: code .
In notebook: Select "Python (qiime)" kernel

Note: qiime environment is separate from the regular .venv environment.

# Notebook 6_Picrust_Functional 

# Activate qiime
conda activate qiime2-2023.7
# install in terminal
pip install --no-deps https://github.com/picrust/q2-picrust2/archive/refs/heads/master.zip
# Refresh the QIIME2 cache:
qiime dev refresh-cache
# Verify the installation:
qiime picrust2 --help
should see the available PICRUSt2 commands listed5.
# If that doesnt work 
conda install -c bioconda -c conda-forge picrust2
pip install --no-deps https://github.com/picrust/q2-picrust2/archive/refs/heads/master.zip
# after installation refresh the QIIME2 cache:
qiime dev refresh-cache
# Verify the installation:
qiime picrust2 --help
# or you can do like me and avoid the drama of the installations by using colab

