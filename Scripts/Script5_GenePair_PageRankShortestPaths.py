#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#FILE NAME: AllGenePairs_PageRank_CentralityMeasures_SupplementaryFilePreparation.py



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 20:46:14 2024

@author: yavuzb2
"""
 
import os
import pandas as pd
import networkx as nx
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
  
def read_gene_pair_pathways(data_folder ='Data', file_name = 'GenePair_ShortestPath_PassThroughOneORMoreThanOneSignalingPathway.txt'):
    """
    Reads a tab-delimited text file containing gene pairs and pathways, and 
    returns three dictionaries based on the number of pathways associated with 
    each gene pair.

    Dictionaries:
        - GenePairDict_OnePathways: Gene pairs with shortest path genes are enriched in exactly one pathway.
        - GenePairDict_TwoPathways: Gene pairs with shortest path genes are enriched in exactly two pathways.
        - GenePairDict_MorethanTwoPathways: Gene pairs with shortest path genes are enriched in more than two pathways.


    Returns:
        tuple: A tuple of three dictionaries (GenePairDict_OnePathways, 
               GenePairDict_TwoPathways, GenePairDict_MorethanTwoPathways).
    """
    # Define the path to the file
    file_path = os.path.join(data_folder, file_name)

    # Initialize dictionaries to store the data based on pathway count
    GenePairDict_OnePathways = {}
    GenePairDict_TwoPathways = {}
    GenePairDict_MorethanTwoPathways = {}

    # Open the file and read it line by line
    with open(file_path, 'r') as file:
        for line in file:
            # Split the line into columns based on tab delimiter
            columns = line.strip().split('\t')
            
            # Extract gene1, gene2, and pathways
            gene1 = columns[0]
            gene2 = columns[1]
            pathways = set(columns[2].split(','))  # Split pathways by comma and store them in a set
            
            # Create a tuple for the gene pair
            gene_pair = (gene1, gene2)
            
            # Store the gene pair in the appropriate dictionary based on pathway count
            pathway_count = len(pathways)
            if pathway_count == 1:
                GenePairDict_OnePathways[gene_pair] = pathways
            elif pathway_count == 2:
                GenePairDict_TwoPathways[gene_pair] = pathways
            else:
                GenePairDict_MorethanTwoPathways[gene_pair] = pathways

    return (GenePairDict_OnePathways, GenePairDict_TwoPathways, GenePairDict_MorethanTwoPathways)



GenePairDict_OnePathways = read_gene_pair_pathways()[0]
GenePairDict_TwoPathways = read_gene_pair_pathways()[1]
GenePairDict_MorethanTwoPathways = read_gene_pair_pathways()[2]


##########
import os
import pandas as pd

def gene_couple_shortest_path_genes(data_folder='Data', file_name='EnrichRAllDiffGenePairs_CouplesShortestPathsSignalingPathway_Info.txt'):
    """
    Reads a tab-delimited text file containing gene couples, signaling pathways, 
    and matching genes, and returns a dictionary where each key is a tuple of 
    (gene1, gene2, pathway) and the value is a set of proteins involved in the shortest path.

    Returns:
        dict: A dictionary with (gene1, gene2, pathway) as keys and sets of 
              shortest path proteins as values.
    """
    # Define the path to the file
    file_path = os.path.join(data_folder, file_name)

    # Read the file into a pandas DataFrame
    df_enrich = pd.read_csv(file_path, sep="\t")

    # Initialize an empty dictionary to store the data
    gene_couple_pathway_shortest_path_proteins = {}

    # Iterate over the DataFrame rows
    for _, row in df_enrich.iterrows():
        # Extract gene1, gene2, pathway, and matching genes
        gene1 = row['Gene1']
        gene2 = row['Gene2']
        pathway = row['SignalingPathway']
        pathway_genes = set(row['MatchingGenesInPathway'].split(";"))

        # Store the gene couple, pathway, and matching genes in the dictionary
        gene_couple_pathway_shortest_path_proteins[(gene1, gene2, pathway)] = pathway_genes

    return gene_couple_pathway_shortest_path_proteins


gene_couple_pathway_shortest_path_proteins = gene_couple_shortest_path_genes()
########

def seed_genes_list(gene_pair, gene_pair_dict_pathway_info):
    """
    Returns the set of seed genes for a given gene pair, including both the 
    source and target nodes, and the proteins involved in the shortest paths 
    for all associated pathways.

    Parameters:
        gene_pair (tuple): A tuple containing the two genes of interest.
        gene_pair_dict_pathway_info (dict): A dictionary with gene pairs as keys 
                                            and associated pathways as values.

    Returns:
        set: A set of seed genes including the source, target, and proteins 
             involved in the shortest path for the given pathways.
    """
    # Get the dictionary of shortest path proteins for each gene couple and pathway
    gene_couple_pathway_shortest_path_proteins = gene_couple_shortest_path_genes()

    # Extract the two genes from the gene pair
    gene1, gene2 = gene_pair

    # Initialize the seed genes set with the source and target genes
    seed_genes = {gene1, gene2}

    # Retrieve pathways associated with the gene pair
    pathways = gene_pair_dict_pathway_info[gene_pair]

    # Ensure pathways is a set to handle both single and multiple pathways uniformly
    if isinstance(pathways, str):
        pathways = {pathways}
    
    # Iterate over each pathway and collect the shortest path proteins
    for pathway in pathways:
        seed_genes.update(gene_couple_pathway_shortest_path_proteins.get((gene1, gene2, pathway), set()))

    return seed_genes

# Example usage
gene_pair = ('MED12', 'PIK3CA')
seed_genes = seed_genes_list(gene_pair, GenePairDict_OnePathways)

# Print the seed genes set (optional)
print(seed_genes)

def create_folder_in_data(folder_name, data_folder='Data'):
    """
    Creates a folder within the specified data directory if it does not already exist.
    """
    # Construct the full path to the folder within the data directory
    folder_path = os.path.join(data_folder, folder_name)
    
    # Create the folder if it does not already exist
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        
create_folder_in_data(folder_name = 'GenePairSpecificSubnetworks-CentralityMeasures')
 
def gene_pair_specific_subnetwork(gene_pair, gene_pair_dict_pathway_info, 
                                  network_file="Data/HIPPI_PPI_network_self_loops_removed.txt", 
                                  output_path='Data/GenePairSpecificSubnetworks-CentralityMeasures'):
    """
    Constructs a gene pair-specific subnetwork based on shortest path enriched pathways 
    and exports it as a Cytoscape SIF file.

    Parameters:
        gene_pair (tuple): A tuple containing the two genes of interest.
        gene_pair_dict_pathway_info (dict): A dictionary with gene pairs as keys and 
                                            associated pathways as values.
        network_file (str): The file name of the PPI network to be used.
        output_path (str): The directory path where the subnetwork SIF file will be saved.
    """
    # Create output directory if it doesn't exist
   

    # Get the seed genes associated with the gene pair
    seed_genes = seed_genes_list(gene_pair, gene_pair_dict_pathway_info)

    # Load the PPI network data
    ppi_data = pd.read_csv(network_file, sep='\t', lineterminator='\n', header=None)
    ppi_data.columns = ['source', 'target']

    # Create a directed graph from the PPI data
    G = nx.from_pandas_edgelist(ppi_data, source='source', target='target', create_using=nx.DiGraph())

    # Initialize personalized PageRank with seed genes
    personalization = {gene: 1 for gene in seed_genes}

    # Calculate PageRank scores
    pagerank_scores = nx.pagerank(G, alpha=0.85, personalization=personalization)

    # Determine threshold for subnetwork construction
    threshold = min(pagerank_scores[gene_pair[0]], pagerank_scores[gene_pair[1]])

    # Select genes for the subnetwork based on the threshold
    if len([gene for gene, score in pagerank_scores.items() if score >= threshold]) <= 50:
        subnetwork_genes = [gene for gene, score in pagerank_scores.items() if score >= threshold - (threshold / 1.2)]
    else:
        subnetwork_genes = [gene for gene, score in pagerank_scores.items() if score >= threshold]

    # Create the subnetwork
    subnetwork = G.subgraph(subnetwork_genes)

    # Export the subnetwork as a Cytoscape SIF file
    added_edges = set()
    sif_filename = os.path.join(output_path, 
                                f'GenePair_{gene_pair[0]}_{gene_pair[1]}_ShortestPathEnrichedPathways_Subnetwork_HIPPIE_ForCentralityMeasures.txt')
    with open(sif_filename, 'w') as sif_file:
        sif_file.write("Node1\tNode2\n")
        for edge in subnetwork.edges():
            if edge[0] != edge[1] and (edge not in added_edges) and ((edge[1], edge[0]) not in added_edges):
                sif_file.write(f"{edge[0]}\t{edge[1]}\n")
                added_edges.add(edge)


Merged_GenePairDict_OneTwo_MoreThanTwoPathways = {}
Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_OnePathways)
Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_TwoPathways)
Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_MorethanTwoPathways)

Merged_GenePairDict_OneTwo_MoreThanTwoPathways[('MED12', 'PIK3CA')]

for gene_pair in Merged_GenePairDict_OneTwo_MoreThanTwoPathways.keys():
    gene_pair_specific_subnetwork(gene_pair,gene_pair_dict_pathway_info = Merged_GenePairDict_OneTwo_MoreThanTwoPathways, network_file =  "Data/HIPPI_PPI_network_self_loops_removed.txt")   












import os
def create_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
create_folder(folder_name = 'GenePairSpecificSubnetworks-CeentralityMeasures')        
       
def GenePairSpecificSubnetwork(gene_pair,GenePairDict_PathwayInformation , network_file =  "Data/HIPPI_PPI_network_self_loops_removed.txt",PATH = 'GenePairSpecificSubnetworks-CeentralityMeasures'):
    ######'HIPPIE-HUMAN-FROM-CX-FILE-NDEX.csv'
    SeedGenes  = seed_genes_list(gene_pair,GenePairDict_PathwayInformation )
    import pandas as pd
   # Import/Create the network that PathLinker will run on
    # create a new network by importing the data from a sample using pandas
    df1 = pd.read_csv(network_file, sep=',', lineterminator='\n')
    df1.columns
   
    df1 = pd.read_csv(network_file, sep='\t', lineterminator='\n',header=None)
    df1.columns = ['source', 'target']
    df= df1
    ppi_data = df [["source","target"]]
    ppi_data.columns =["source","target"]  
    pair = gene_pair
    Seed_Dict ={k:1 for k in SeedGenes.difference({pair[0],pair[1]})}    
    Seed_Dict[pair[0]] = 2
    Seed_Dict[pair[1]] = 2
    #TissueSpecificSeedGenes is the new set of seed genes
    import pandas as pd
    import networkx as nx
    # Create a directed graph from the PPI data
    G = nx.from_pandas_edgelist(ppi_data, source='source', target='target', create_using=nx.DiGraph())

    # List of seed genes
    seed_genes = SeedGenes  # Replace with  seed genes

    # Initialize personalized PageRank with seed genes
    personalization = {gene: 1 for gene in seed_genes}

    # Calculate PageRank
    pagerank_scores = nx.pagerank(G, alpha=0.85, personalization=personalization) #alpha=0.85,
    

    # Set a threshold for subnetwork construction
    #threshold = 0.002  # Adjust the threshold as needed, breast0.002
    threshold = min(pagerank_scores[pair[0]], pagerank_scores[pair[1]]) #BREAST 0.002  set threshold as minimum(pagerank_scores["ESR1"], pagerank_scores["PIK3CA"])
    # Create a subnetwork by including genes with PageRank scores above the threshold
    if len([gene for gene, score in pagerank_scores.items() if score >= threshold])<=50:
        subnetwork_genes = [gene for gene, score in pagerank_scores.items() if score >= threshold-(threshold/1.2)] #to reduce the numbers of subnetwork nodes
    else:
        subnetwork_genes = [gene for gene, score in pagerank_scores.items() if score >= threshold]
    subnetwork = G.subgraph(subnetwork_genes)


    # Export su8bnetwork as Cytoscape SIF file
    ADD = set()
    sif_filename = PATH+'GenePair_{}_{}_ShortestPathEnrichedPathways_Subnetwork_HIPPIE_ForCentralityMEasures.txt'.format(pair[0],pair[1])
    with open(sif_filename, 'a') as sif_file:
        sif_file.write("Node1" +"\t"+"Node2"+"\n")
    with open(sif_filename, 'a') as sif_file:
        for edge in subnetwork.edges():
            if edge[0] != edge[1]:
                    if ((edge[0],edge[1]) not in ADD) and ((edge[1],edge[0]) not in ADD):
                        sif_file.write(str(edge[0])+"\t"+str(edge[1])+"\n")
                        ADD.add((edge[0],edge[1]))
                        ADD.add((edge[1],edge[0]))
    ADD           


#Call the above function

#GenePairSpecificSubnetwork(gene_pair = ("ESR1","PIK3CA"),GenePairDict_PathwayInformation =GenePairDict_TwoPathways, network_file =  "HIPPI_PPI_network_self_loops_removed.txt")   



Merged_GenePairDict_OneTwo_MoreThanTwoPathways = {}
Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_OnePathways)
Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_TwoPathways)
Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_MorethanTwoPathways)

Merged_GenePairDict_OneTwo_MoreThanTwoPathways[('MED12', 'PIK3CA')]

for gene_pair in Merged_GenePairDict_OneTwo_MoreThanTwoPathways.keys():
    GenePairSpecificSubnetwork(gene_pair,GenePairDict_PathwayInformation = Merged_GenePairDict_OneTwo_MoreThanTwoPathways, network_file =  "Data/HIPPI_PPI_network_self_loops_removed.txt")   


# import networkx as nx
# def calculate_centrality_measures(subnetwork):
    
#     # Calculate degree centrality
#     degree_centrality = nx.degree_centrality(subnetwork)
    
#     # Calculate betweenness centrality
#     betweenness_centrality = nx.betweenness_centrality(subnetwork)
    
#     # Calculate closeness centrality
#     closeness_centrality = nx.closeness_centrality(subnetwork)
    
#     # Calculate eigenvector centrality
#     eigenvector_centrality = nx.eigenvector_centrality(subnetwork)
    
#     return {
#         'degree_centrality': degree_centrality,
#         'betweenness_centrality': betweenness_centrality,
#         'closeness_centrality': closeness_centrality,
#         'eigenvector_centrality': eigenvector_centrality
#     }



# def calculate_Q3(dictionary):
#     """
#     Calculate the interquartile range (IQR) for a list of numbers.

#     Parameters:
#     numbers (list): A list of numeric values.

#     Returns:
#     float: The interquartile range of the numbers.
#     """
    
#     numbers = list(dictionary.values())
#     sorted_numbers = sorted(numbers)
#     n = len(sorted_numbers)

#     def percentile(data, percent):
#         k = (len(data) - 1) * percent / 100.0
#         f = int(k)
#         c = k - f
#         if f + 1 < len(data):
#             return data[f] + (data[f + 1] - data[f]) * c
#         else:
#             return data[f]

#     q1 = percentile(sorted_numbers, 25)
#     q3 = percentile(sorted_numbers, 75)
    
#     return  q3


 

def calculate_centrality_measures(subnetwork):
    
    # Calculate degree centrality
    degree_centrality = nx.degree_centrality(subnetwork)
    
    # Calculate betweenness centrality
    betweenness_centrality = nx.betweenness_centrality(subnetwork)
    
    # Calculate closeness centrality
    closeness_centrality = nx.closeness_centrality(subnetwork)
    
    # Calculate eigenvector centrality
    eigenvector_centrality = nx.eigenvector_centrality(subnetwork)
    
    return {
        'degree_centrality': degree_centrality,
        'betweenness_centrality': betweenness_centrality,
        'closeness_centrality': closeness_centrality,
        'eigenvector_centrality': eigenvector_centrality
    }



def calculate_Q3(dictionary):
    """
    Calculate the interquartile range (IQR) for a list of numbers.

    Parameters:
    numbers (list): A list of numeric values.

    Returns:
    float: The interquartile range of the numbers.
    """
    
    numbers = list(dictionary.values())
    sorted_numbers = sorted(numbers)
    n = len(sorted_numbers)

    def percentile(data, percent):
        k = (len(data) - 1) * percent / 100.0
        f = int(k)
        c = k - f
        if f + 1 < len(data):
            return data[f] + (data[f + 1] - data[f]) * c
        else:
            return data[f]

    q1 = percentile(sorted_numbers, 25)
    q3 = percentile(sorted_numbers, 75)
    
    return  q3


def SelectTargetsWithBetweennessCentrality(gene_pair,PATH = 'GenePairSpecificSubnetworks-CeentralityMeasures'):

    prefixed = PATH+'GenePair_{}_{}_ShortestPathEnrichedPathways_Subnetwork_HIPPIE_ForCentralityMEasures.txt'.format(gene_pair[0],gene_pair[1])
    Subnetwork_Edges = []
    with open(prefixed,"r") as infile:
        for line in infile:
            if ('Node' not in line) and ('nan' not in line):
                splitted = line.rstrip("\n").split("\t")
                Subnetwork_Edges.append((splitted[0],splitted[1]))
            
  
    # Create an example subnetwork
    subnetwork = nx.Graph()
    
    # Adding nodes and edges (for demonstration purposes, replace with your actual data)
    subnetwork.add_edges_from(Subnetwork_Edges)
    
    # Calculate centrality measures
    centrality_measures = calculate_centrality_measures(subnetwork)
    # Example usage
    my_dict = centrality_measures['betweenness_centrality']

    q3 = calculate_Q3(my_dict)
    threshold = q3
    
    
    DrugTargets = set()
    for gene in my_dict.keys():
        if my_dict[gene] > threshold:
            DrugTargets.add(gene)
            
    return DrugTargets

for pair in Merged_GenePairDict_OneTwo_MoreThanTwoPathways.keys():
    Targets = SelectTargetsWithBetweennessCentrality(pair,PATH = 'Data/GenePairSpecificSubnetworks-CentralityMeasures/')
    target = '+'.join(list(Targets))
    with open('Data/GenePair_DrugTargetsFromBetweennessCentrality_IQR_GreaterThan_Q3.txt','a') as outfile:
        outfile.write(pair[0]+"\t"+pair[1]+"\t"+target+'\n')
   
#######################   

#convert this to function below
    
# GenePairsTwo  = [pair[0]+"|"+pair[1] for pair in GenePairDict_TwoPathways.keys()]
# GenePairsThree  = [pair[0]+"|"+pair[1] for pair in GenePairDict_MorethanTwoPathways.keys()]
# GenePairsOne  = [pair[0]+"|"+pair[1] for pair in GenePairDict_OnePathways.keys()]

# DoubletComponentsTwo = set()
# for double in GenePairDict_TwoPathways.keys():
#     DoubletComponentsTwo.add(double[0])
#     DoubletComponentsTwo.add(double[1])
    
# DoubletComponentsThree = set()
# for double in GenePairDict_TwoPathways.keys():
#     DoubletComponentsThree.add(double[0])
#     DoubletComponentsThree.add(double[1])
# DoubletComponentsOne = set()
# for double in GenePairDict_OnePathways.keys():
#     DoubletComponentsOne.add(double[0])
#     DoubletComponentsOne.add(double[1])
        
# ###########
# #############################

# # Step 1: Read data from a text file into a pandas DataFrame
# filename = 'Data/GenePair_DrugTargetsFromBetweennessCentrality_IQR_GreaterThan_Q3.txt'#


# # Step 1: Read data from a text file into a pandas DataFrame

# # Assuming your data file is tab-separated with columns GeneA, GeneB, Interactions
# df1 = pd.read_csv(filename, sep='\t', header=None, names=['GeneA', 'GeneB', 'Interactions'])

# # Step 2: Merge GeneA and GeneB into a single 'GenePair' column
# df1['GenePair'] = df1['GeneA'] + '|' + df1['GeneB']

# df2_Two = df1[df1['GenePair'].isin(GenePairsTwo)]
# df2_Three = df1[df1['GenePair'].isin(GenePairsThree)]
# df2_One = df1[df1['GenePair'].isin(GenePairsOne)]


# # Drop columns with all 0 values and columns where the sum is less than 5

# # Step 3: Create a set of all unique interactions
# #175 drug target genes
# interactionsTwo = set()
# for interaction_list in df2_Two['Interactions']:
#     interactionsTwo.update(interaction_list.split('+'))
# interactionsTwo = interactionsTwo.union(DoubletComponentsTwo)

# #175 drug target genes
# interactionsThree = set()
# for interaction_list in df2_Three['Interactions']:
#     interactionsThree.update(interaction_list.split('+'))
# interactionsThree = interactionsThree.union(DoubletComponentsThree)

# #175 drug target genes
# interactionsOne = set()
# for interaction_list in df2_One['Interactions']:
#     interactionsOne.update(interaction_list.split('+'))
# interactionsOne = interactionsOne.union(DoubletComponentsOne)

# intersection_of_all_connector_nodes = interactionsTwo.intersection(interactionsThree,interactionsOne)  
   
###update function

import pandas as pd

def process_gene_pair_interactions(gene_pair_dict_one, gene_pair_dict_two, gene_pair_dict_three, filename):
    """
    Processes gene pair interactions by merging gene pairs, filtering data based on gene pairs,
    and calculating the intersection of all connector nodes.

    Parameters:
        gene_pair_dict_one (dict): Dictionary with gene pairs and associated pathways (one pathway).
        gene_pair_dict_two (dict): Dictionary with gene pairs and associated pathways (two pathways).
        gene_pair_dict_three (dict): Dictionary with gene pairs and associated pathways (more than two pathways).
        filename (str): Path to the input file containing gene pair interactions.

    Returns:
        set: Intersection of all connector nodes across the three gene pair sets.
    """
    # Helper function to merge gene pairs into a single string
    def merge_gene_pairs(gene_pair_dict):
        return [f"{pair[0]}|{pair[1]}" for pair in gene_pair_dict.keys()]
    
    # Generate lists of merged gene pairs
    gene_pairs_one = merge_gene_pairs(gene_pair_dict_one)
    gene_pairs_two = merge_gene_pairs(gene_pair_dict_two)
    gene_pairs_three = merge_gene_pairs(gene_pair_dict_three)

    # Helper function to extract unique gene components from gene pairs
    def extract_gene_components(gene_pair_dict):
        components = set()
        for gene1, gene2 in gene_pair_dict.keys():
            components.add(gene1)
            components.add(gene2)
        return components

    # Extract unique gene components for each pathway category
    components_one = extract_gene_components(gene_pair_dict_one)
    components_two = extract_gene_components(gene_pair_dict_two)
    components_three = extract_gene_components(gene_pair_dict_three)

    # Read data from the input file
    df = pd.read_csv(filename, sep='\t', header=None, names=['GeneA', 'GeneB', 'Interactions'])

    # Merge GeneA and GeneB into a single 'GenePair' column
    df['GenePair'] = df['GeneA'] + '|' + df['GeneB']

    # Filter the DataFrame based on the gene pairs
    df_filtered_one = df[df['GenePair'].isin(gene_pairs_one)]
    df_filtered_two = df[df['GenePair'].isin(gene_pairs_two)]
    df_filtered_three = df[df['GenePair'].isin(gene_pairs_three)]

    # Helper function to extract interactions from a DataFrame and merge with gene components
    def extract_interactions(df_filtered, components):
        interactions = set()
        for interaction_list in df_filtered['Interactions']:
            interactions.update(interaction_list.split('+'))
        return interactions.union(components)

    # Extract interactions for each pathway category and merge with components
    interactions_one = extract_interactions(df_filtered_one, components_one)
    interactions_two = extract_interactions(df_filtered_two, components_two)
    interactions_three = extract_interactions(df_filtered_three, components_three)

    # Calculate the intersection of all connector nodes across the three sets
    intersection_of_all_connector_nodes = interactions_one.intersection(interactions_two, interactions_three)

    return intersection_of_all_connector_nodes

# Example usage
filename = 'Data/GenePair_DrugTargetsFromBetweennessCentrality_IQR_GreaterThan_Q3.txt'
intersection_nodes = process_gene_pair_interactions(GenePairDict_OnePathways, 
                                                   GenePairDict_TwoPathways, 
                                                   GenePairDict_MorethanTwoPathways, 
                                                   filename)

print(intersection_nodes)
  
len(intersection_nodes)
