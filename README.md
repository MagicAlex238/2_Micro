# MIC_Corrosion
This is the repository for the Dissertation "Nobel Microorganisms associated with MIC"
The first notebook is a visualisation the**Physicochemical** data analysis. These repository corresponds to the Session 3.3 and 4.5 in the main text.

| Notebook                    |Description                                                                                                  |
|----------------------------|--------------------------------------------------------------------------------------------------------------|
|**1_Splitdf.ipynb**     | Explaination of the data, then I separate the data into categories according to the failure study and the physicochemical study that distinguise three differnet degrees of corrosion.  The notebook analyse the type of data distribution. More details can be found in the main text of the dissertation (4.4.1). |
| **2_Filtering.ipynb**     | A notebook 1. Identify the most relevant bacteria genus(GIDs) amongts the 880 genus by filering out the percentage that is in a concentration of less than 2%, the genera are then sorted and an extracolumn tells the influence frecuency of those genera, lastly the data is divided into three traffic lights. A small statistic is done with the number of GID, count, Unique, top and frequency. |
| **3_Predominance.ipynb**      | An example script to determine predominant chemical species. Refer to the main text (methods 3.3.5) and the results and discussion (Chapter 4.4.6) for details. |
| **4_Sensitivity.ipynb**      | A notebook for sensitivity analysis to assess the influence of organic matter parameters on chemical species distribution. |



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
