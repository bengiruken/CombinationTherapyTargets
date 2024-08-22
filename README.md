## Identifying Target Proteins for Combination Therapies 

CombinationTherapyTargets is an innovative network-based approach designed to identify proteins suitable for combination cancer therapies. 

To temper the expected resistance to single drug regimen, we offer a concept-based stratified pipeline aimed at selecting co-targets for drug combinations. Our strategy is unique in its co-target selection being based on signaling pathways. This is significant since in cancer, drug resistance commonly bypasses blocked proteins by wielding alternative, or complementary, routes to execute cell proliferation. Our network-informed signaling-based approach harnesses advanced network concepts and metrics, and our compiled tissue-specific co-existing mutations. Co-existing driver mutations are common in resistance. Thus, to mimic cancer and counter drug resistance scenarios, our pipeline seeks co-targets that when targeted by drug combinations, can shut off cancer’s modus operandi. That is, its parallel or complementary signaling pathways would be blocked. Rotating through combinations could further lessen emerging resistance. 


## Installation

git clone https://github.com/bengiruken/CombinationTherapyTargets.git
cd CombinationTherapyTargets
pip install .

## Usage
To run the tool and execute all the modules, follow these steps:

If you haven’t already, ensure you’re in the project directory:
cd CombinationTherapyTargets
Use Python to run the main.py file:
python main.py

## Requirements
Requirements

The following Python packages are required and will be installed automatically:

pandas==2.0.3
scipy==1.11.1
numpy==1.24.3
gseapy==1.1.2
matplotlib-inline==0.1.6
seaborn==0.12.2
networkx==3.2.1



## Data

### List of gene pairs harboring co-existing mutations

In this study, previously discovered co-exisiting mutations were retrieved from tumor-specific networks were retrieved from [(Yavuz et al. 2024)](doi: https://doi.org/10.1101/2024.05.01.592039). 

### HIPPIE PPI network

HIPPE protein-protein interaction (PPI) network is downloaded from (https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/).


### Shortest paths
Shortest paths connecting gene pairs with co-exisiting mutations were computed by using PathLinker algorithm (https://github.com/Murali-group/PathLinker). The corresponding output is in the paths.txt file.


### KEGG pathways

KEGG_2019_human signaling pathways were obtained from EnrichR library (https://maayanlab.cloud/Enrichr/#libraries)


### Transcription factors 
The list of transcription factors were obtained from TRRUST (https://www.grnpedia.org/trrust/).

### Oncogenes and tumor suppressor genes
The list of oncogenes and tumor suppressor genes were obtained from OncoKB (https://www.oncokb.org/cancer-genes).

## Scripts

All the corresponding analysis conducted with Python language and the corresponding data were stored in [this GitHub repository](https://github.com/bengiruken/CombinationTherapyTargets).

### Structure

There are two folders having scripts and data filesnamed as **"Scripts"** and **"Data"**, respectively.

CombinationTherapyTargets file should be downloaded and main.py module needs to be compiled. 


