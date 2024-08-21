#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:32:16 2024

@author: yavuzb2
"""

'AfterBiowulf_EnrichR_Result_Processing.py'

#import functions for data curation
from Scripts import Script1_FunctionsForDataCuration  as DC

list_of_all_gene_pairs = DC.ListOfGenePairsWithCoexisitngMutations()


import pandas as pd
import os

def get_cooccurring_gene_couples(file_path='Data', file_name="All_Different_Gene_Doublets_Supplementary_Table1.xlsx", tendency='Co-occurence', min_double_mut=5):
    """
    Get the list of co-occurring gene doublets from an Excel file.

    Args:
        file_path (str): Path to the directory containing the Excel file.
        file_name (str): Name of the Excel file.
        tendency (str): The tendency filter for selecting gene pairs. Default is 'Co-occurence'.
        min_double_mut (int): Minimum number of double mutations required to consider a gene pair. Default is 5.

    Returns:
        set: A set of tuples representing co-occurring gene pairs.
    """
    # Construct the full path to the Excel file
    file_full_path = os.path.join(file_path, file_name)
    
    # Read the Excel file into a DataFrame
    df_met = pd.read_excel(file_full_path)
    
    # Filter the DataFrame based on the specified tendency and minimum double mutations
    df_met_filtered = df_met[(df_met["Tendency"] == tendency) & (df_met['DoubleMut#'] >= min_double_mut)]
    df_met_filtered.reset_index(drop=True, inplace=True)
    
    # Extract the list of co-occurring gene pairs
    cooccurring_gene_couples = set()
    for ind in df_met_filtered.index:
        gene1 = df_met_filtered.loc[ind, "Mut1"].split("_")[0]
        gene2 = df_met_filtered.loc[ind, "Mut2"].split("_")[0]
        cooccurring_gene_couples.add((gene1, gene2))
    
    return cooccurring_gene_couples

if __name__ == "__main__":
    # Example usage of the function
    gene_couples = get_cooccurring_gene_couples()
    print(f"Co-occurring gene couples: {gene_couples}")


import pandas as pd
import os
import matplotlib.pyplot as plt

file_path='Data'
file_name='EnrichRAllDiffGenePairs_CouplesShortestPathsSignalingPathway_Info.txt'

# Construct the full path to the input file
file_full_path = os.path.join(file_path, file_name)

# Read the EnrichR results into a DataFrame
df_enrich = pd.read_csv(file_full_path, sep="\t")



# Split the 'MatchInPathway' column into 'Numerator' and 'Denominator'
df_enrich[['Numerator', 'Denominator']] = df_enrich['MatchInPathway'].str.split('/', expand=True)

print(df_enrich[['Numerator', 'Denominator']])

# Convert columns to numeric and handle errors
df_enrich['Numerator'] = pd.to_numeric(df_enrich['Numerator'], errors='coerce')
df_enrich['Denominator'] = pd.to_numeric(df_enrich['Denominator'], errors='coerce')
df_enrich['SizeshortesPath'] = pd.to_numeric(df_enrich['SizeshortesPath'], errors='coerce')

# Print to check the converted columns
print(df_enrich[['Numerator', 'Denominator','SizeshortesPath']])
# Convert to numeric and coerce errors to NaN
# df_enrich['Numerator'] = pd.to_numeric(df_enrich['Numerator'], errors='coerce')
# df_enrich['Denominator'] = pd.to_numeric(df_enrich['Denominator'], errors='coerce')

df_enrich['Numerator'] = df_enrich['Numerator'].astype(float)
df_enrich['Denominator'] = df_enrich['Denominator'].astype(float)
 
 # Calculate the fraction
df_enrich['Fraction'] = (df_enrich['Numerator'] / df_enrich['SizeshortesPath']) * 100
 

def process_enrichr_results(file_path='Data', file_name='EnrichRAllDiffGenePairs_CouplesShortestPathsSignalingPathway_Info.txt', threshold=20):
    """
    Process EnrichR results to extract and analyze gene pairs and their signaling pathways.

    Args:
        file_path (str): Path to the directory containing the input file.
        file_name (str): Name of the input file.
        threshold (float): Threshold for the fraction to consider a pathway significant.

    Returns:
        dict: A dictionary where keys are gene pairs and values are sets of significant signaling pathways.
    """
    # Construct the full path to the input file
    file_full_path = os.path.join(file_path, file_name)
    
    # Read the EnrichR results into a DataFrame
    df_enrich = pd.read_csv(file_full_path, sep="\t")
    
    
    
    # Split the 'MatchInPathway' column into 'Numerator' and 'Denominator'
    df_enrich[['Numerator', 'Denominator']] = df_enrich['MatchInPathway'].str.split('/', expand=True)
    
    # # Convert to numeric and coerce errors to NaN
    # df_enrich['Numerator'] = pd.to_numeric(df_enrich['Numerator'], errors='coerce')
    # df_enrich['Denominator'] = pd.to_numeric(df_enrich['Denominator'], errors='coerce')
    # Convert columns to numeric and handle errors
    df_enrich['Numerator'] = pd.to_numeric(df_enrich['Numerator'], errors='coerce')
    df_enrich['Denominator'] = pd.to_numeric(df_enrich['Denominator'], errors='coerce')
    df_enrich['SizeshortesPath'] = pd.to_numeric(df_enrich['SizeshortesPath'], errors='coerce')

    # Print to check the converted columns
    print(df_enrich[['Numerator', 'Denominator','SizeshortesPath']])
    # Convert to numeric and coerce errors to NaN
    # df_enrich['Numerator'] = pd.to_numeric(df_enrich['Numerator'], errors='coerce')
    # df_enrich['Denominator'] = pd.to_numeric(df_enrich['Denominator'], errors='coerce')

    df_enrich['Numerator'] = df_enrich['Numerator'].astype(float)
    df_enrich['Denominator'] = df_enrich['Denominator'].astype(float)
     
     # Calculate the fraction
    df_enrich['Fraction'] = (df_enrich['Numerator'] / df_enrich['SizeshortesPath']) * 100
    
    # Plot the histogram of the 'Fraction' column
    ax = df_enrich['Fraction'].plot.hist(bins=12, alpha=0.5)
    plt.xlabel('Fraction (%)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Fractions')
    plt.show()
    
    # Initialize a dictionary to store gene pairs and their significant pathways
    gene_pair_signal_passing_pathways = {}
    
    # Iterate through the DataFrame to populate the dictionary
    for ind in df_enrich.index:
        gene1 = df_enrich.loc[ind, 'Gene1']
        gene2 = df_enrich.loc[ind, 'Gene2']
        frac = df_enrich.loc[ind, 'Fraction']
        pathway = df_enrich.loc[ind, 'SignalingPathway']
        gene_pair = (gene1, gene2)
        
        if gene_pair not in gene_pair_signal_passing_pathways:
            gene_pair_signal_passing_pathways[gene_pair] = set()
        
        if frac > threshold:
            gene_pair_signal_passing_pathways[gene_pair].add(pathway)
    
    return  gene_pair_signal_passing_pathways#df_enrich

if __name__ == "__main__":
    gene_pair_pathways = process_enrichr_results()
    print(gene_pair_pathways.get(("ESR1", "PIK3CA")))
    print(gene_pair_pathways.get(("PIK3CA", "PTEN")))

# #there are 1263 gene pairs with shortest path information
# Zero =set()
# for pair in gene_pair_pathways.keys():
#     pathways = gene_pair_pathways[pair]
#     if len(pathways) == 0 :
#         Zero.add(pair)
#1091 gene pair's shortest paths are not enriched in a signaling pathway

''' ***********  '''



def write_signaling_pathways_to_file( output_path='Data'):
    """
    Write gene pairs and their signaling pathways to files based on the number of pathways they are enriched in.

    Args:
        gene_pair_signal_passing_pathways (dict): Dictionary of gene pairs and their signaling pathways.
        result_matrix (pd.DataFrame): DataFrame containing similarity scores between pathways.
        output_path (str): Path to the directory where output files will be saved.
    """
    signaling_file = os.path.join(output_path, 'GenePair_ShortestPath_PassThroughOneORMoreThanOneSignalingPathway.txt')
    gene_pair_pathways = process_enrichr_results()
    for pair in gene_pair_pathways.keys():
        pathways = gene_pair_pathways[pair]
        if len(pathways) > 0 :
            NumberOfPathways = str(len(pathways))
            PATHWAYS = ",".join(list(pathways))
            with open(signaling_file, "a") as outfile:
                outfile.write(f"{pair[0]}\t{pair[1]}\t{PATHWAYS}\t{NumberOfPathways}\n")
                                   


if __name__ == "__main__":
    gene_pair_pathways = process_enrichr_results()
    write_signaling_pathways_to_file()




















