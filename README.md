## Discovering Anticancer Drug Target Combinations via Network-Informed Signaling-Based Approach 

Cancer requires drug treatment. The oncologist is faced with the problem: which drugs to prescribe and which regimen? Much is still unknown, and the number of possible drugs, and their combinations is vast. At the outset, the oncologist reckons with at least two established facts: (i) patients receiving successive single molecules treatments are likely to experience drug resistance, and (ii), to select optimal drug combinations requires to pick the ‘best’ protein drug target combinations. Intuitively, target selection should precede drug selection, implying that well-informed strategies would opt to first consider drug targets — not drugs — combinations. Nowadays, drug combinations that oncologists consider are empirical and limited. They are restricted primarily by observations and praxis, that is, scant clinical experience with their application. 
Here we develop a strategy for selecting optimal drug target combinations following nature. We use protein-protein interaction networks and shortest paths to discover communication pathways in cells based on interaction network topology. Our strategy mimics cancer signaling in drug resistance, which commonly harnesses pathways parallel to those blocked by drugs, thereby bypassing them. 
We select key communication nodes as combination drug targets inferred from topological features of networks. We test our network-informed signaling-based approach to discover anticancer drug target combinations on available clinical data, patient-derived breast and colorectal cancers. Alpelisib + LJM716 and alpelisib + cetuximab + encorafenib combinations diminish tumors in breast and colorectal cancers, respectively.
Our network-based approach discovers optimal protein co-target combinations to counter resistance, selecting co-targets from alternative pathways and their connectors.


Figure 1 shows the schematic of our rationale:

![FIGURE~1](https://github.com/user-attachments/assets/e83b183a-419b-433e-bcdc-c09fc89c418c)

Schematic of our pipeline:

![FIGURE~3](https://github.com/user-attachments/assets/2ab53205-5512-4552-8d59-e574c74bc34a)

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


