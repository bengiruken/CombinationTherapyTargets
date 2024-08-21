#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:53:44 2024

@author: yavuzb2
"""


"""
Created on Tue Jul 30 14:58:43 2024

@author: yavuzb2

                                COMBINATION THERAPY TARGETS
                                
                                        DATA CURATION
                                        
                                    Bengi Ruken Yavuz, PhD
"""




def ListOfGenePairsWithCoexisitngMutations(file_path='Data', file_name='PanCancer_GenePairs_CoexistingMut_PathLinkerInput.txt'):
    """
    Get the list of gene pairs harboring co-existing mutations.
    
    Parameters:
    file_path (str): The directory path where the file is located.
    file_name (str): The name of the file containing gene pairs with co-existing mutations.

    Returns:
    set: A set of tuples representing gene pairs with co-existing mutations.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    file_full_path = os.path.join(file_path, file_name)
    with open(file_full_path, 'r') as infile:
        lines = infile.readlines()

    # Create a set of tuples
    gene_pairs = {tuple(line.strip().split()) for line in lines}

    return gene_pairs


def HIPPIE_PPI_Network(url, file_path='Data'):
    import os
    import pandas as pd

    file_path = 'Data'
    """
    Download data from the provided URL and curate the HIPPIE PPI network for PathLinker.
    
    Parameters:
    url (str): The URL to download the HIPPIE PPI network data.
    file_path (str): The directory path where the HIPPIE PPI network file is located.
    """
    network_file = 'HIPPIE-HUMAN-FROM-CX-FILE-NDEX.csv'
    file_full_path = os.path.join(file_path, network_file)

    # Create a new network by importing the data from a sample using pandas
    df1 = pd.read_csv(file_full_path, sep=',', lineterminator='\r')
    df1[['source', 'target']] = df1.name.str.split(' (interacts-with) ', expand=True, regex=False)
    df = df1[['source', 'target']]

    # Remove self-loops and write to output file
    output_file = 'HIPPI_PPI_network_self_loops_removed.txt'
    output_full_path = os.path.join(file_path, output_file)
    with open(output_full_path, 'a') as outfile:
        for _, row in df.iterrows():
            source = row['source']
            target = row['target']
            if source != target and pd.notna(source) and pd.notna(target):
                outfile.write(f"{source}\t{target}\n")


def KEGG2019_Pathway_Gene_Dictionary(file_name="KEGG_2019_Human.txt", file_path='Data'):
    """
    Curate the KEGG 2019 signaling pathways.
    
    Parameters:
    file_name (str): The name of the KEGG 2019 dataset file.
    file_path (str): The directory path where the KEGG 2019 dataset file is located.

    Returns:
    dict: A dictionary where keys are pathways and values are sets of genes in those pathways.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    file_full_path = os.path.join(file_path, file_name)
    Pathway_Gene_Dict = {}
    Genes = set()

    with open(file_full_path, "r") as infile:
        for line in infile:
            if "signaling" in line and "pathway" in line:
                splitted = line.rstrip("\n").split("\t")
                genes = set(filter(None, splitted[1:]))
                Pathway_Gene_Dict[splitted[0]] = genes
                Genes.update(genes)
    
    return Pathway_Gene_Dict

def Gene_Pathway_Dictionary( file_name="KEGG_2019_Human.txt", file_path='Data'):
    """
    Create a dictionary mapping genes to pathways.
    
    Parameters:
    pathway_gene_dict (dict): A dictionary where keys are pathways and values are sets of genes.

    Returns:
    dict: A dictionary where keys are genes and values are sets of pathways containing those genes.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    pathway_gene_dict = KEGG2019_Pathway_Gene_Dictionary()
    gene_pathway_dict = {}

    for pathway, genes in pathway_gene_dict.items():
        for gene in genes:
            if gene not in gene_pathway_dict:
                gene_pathway_dict[gene] = set()
            gene_pathway_dict[gene].add(pathway)
    
    return gene_pathway_dict
len(Gene_Pathway_Dictionary()["PIK3CA"])

def TRRUST_Transcription_Factor_Set(file_name="Trrust_rawdata_human.tsv", file_path='Data'):
    """
    Curate the set of transcription factors from TRRUST.

    Parameters:
    file_name (str): The name of the TRRUST dataset file.
    file_path (str): The directory path where the TRRUST dataset file is located.

    Returns:
    set: A set of transcription factors.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    file_full_path = os.path.join(file_path, file_name)
    transcription_factors = set()

    with open(file_full_path, "r") as infile:
        for line in infile:
            splitted = line.split("\t")
            transcription_factors.add(splitted[0])
    
    return transcription_factors


def OncoKB_OG_TSG(file_name="OncoKB_CancerGeneList_OG_TSG_Label.txt", file_path='Data'):
    """
    Curate the set of oncogenes and tumor suppressors from OncoKB.

    Parameters:
    file_name (str): The name of the OncoKB dataset file.
    file_path (str): The directory path where the OncoKB dataset file is located.

    Returns:
    dict: A dictionary with keys 'OG' and 'TSG' containing sets of oncogenes and tumor suppressors respectively.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    oncogenes = set()
    tumor_suppressors = set()
    file_full_path = os.path.join(file_path, file_name)

    with open(file_full_path, "r") as infile:
        for line in infile:
            line = line.rstrip("\n")
            splitted = line.split("\t")
            if splitted[-1] == 'OG':
                oncogenes.add(splitted[0])
            elif splitted[-1] == 'TSG':
                tumor_suppressors.add(splitted[0])
    
    return {'OG': oncogenes, 'TSG': tumor_suppressors}


def RTKs(file_name="ReceptorKinaseList.csv", file_path='Data'):
    """
    Curate the set of Receptor Tyrosine Kinases (RTKs).

    Parameters:
    file_name (str): The name of the RTK dataset file.
    file_path (str): The directory path where the RTK dataset file is located.

    Returns:
    set: A set of RTKs.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    file_full_path = os.path.join(file_path, file_name)
    df = pd.read_csv(file_full_path, sep=",")
    rtks = set(df['Approved symbol'].to_list())
    
    return rtks

import pandas as pd

def CMAP_drug_target_dictionary(drug_data_file='Data/OncologyLaunchedDrugsNetworkInput_CMAP.txt'):
    """
    Creates a dictionary of drug targets from the drug data file.

    Parameters:
    drug_data_file (str): Path to the file containing drug target information.

    Returns:
    dict: A dictionary where keys are drug target proteins and values are lists of drugs targeting those proteins.
    """
    import os
    import pandas as pd

    file_path = 'Data'
    # Read the drug target data into a DataFrame
    drug_targets_df = pd.read_csv(drug_data_file, sep="\t", header=None, names=['Drug', 'Target', 'Mechanism'])

    # Initialize an empty dictionary to store drug targets
    drug_target_dict = {}

    # Iterate through the DataFrame to populate the dictionary
    for _, row in drug_targets_df.iterrows():
        target = row['Target']
        drug = row['Drug']

        if target in drug_target_dict:
            drug_target_dict[target].append(drug)
        else:
            drug_target_dict[target] = [drug]

    return drug_target_dict




# if __name__ == "__main__":
#     # Example usage of the functions
#     gene_pairs = ListOfGenePairsWithCoexisitngMutations()
#     HIPPIE_PPI_Network(url='https://www.ndexbio.org/viewer/networks/89dd3925-3718-11e9-9f06-0ac135e8bacf')
#     kegg_pathway_gene_dict = KEGG2019_Pathway_Gene_Dictionary()
#     transcription_factors = TRRUST_Transcription_Factor_Set()
#     oncogenes_and_tumor_suppressors = OncoKB_OG_TSG()
#     rtks = RTKs()
    
#     # Print results or perform further processing
#     print(f"Gene pairs: {gene_pairs}")
#     print(f"KEGG Pathways: {kegg_pathway_gene_dict}")
#     print(f"Transcription Factors: {transcription_factors}")
#     print(f"Oncogenes and Tumor Suppressors: {oncogenes_and_tumor_suppressors}")
#     print(f"RTKs: {rtks}")
