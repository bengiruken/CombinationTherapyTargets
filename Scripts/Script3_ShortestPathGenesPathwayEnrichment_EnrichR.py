#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 07:54:42 2024

Author: yavuzb2
"""

import os
import gseapy
import numpy as np

def get_gene_couples(file_path, file_name):
    """
    Extracts gene couples with shortest path information from the specified file.

    Parameters:
    file_path (str): The directory path where the file is located.
    file_name (str): The name of the PathLinker output file.

    Returns:
    set: A set of tuples, where each tuple contains a pair of genes.
    """
    gene_pairs = set()
    file_full_path = os.path.join(file_path, file_name)
    
    try:
        with open(file_full_path, 'r') as infile:
            for line in infile:
                splitted = line.rstrip("\n").split("\t")
                path_info = splitted[-1].split("|")
                gene_pair = (path_info[0], path_info[-1])
                gene_pairs.add(gene_pair)
    except FileNotFoundError:
        print(f"Error: File '{file_full_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return gene_pairs


def lump_shortest_path_genes(file_path, file_name, gene_pairs):
    """
    Aggregates all shortest path genes for each gene couple.

    Parameters:
    file_path (str): The directory path where the file is located.
    file_name (str): The name of the PathLinker output file.
    gene_pairs (set): A set of gene pairs.

    Returns:
    dict: A dictionary where keys are gene pairs and values are sets of shortest path genes.
    """
    couple_shortest_path_dict = {couple: set() for couple in gene_pairs}
    file_full_path = os.path.join(file_path, file_name)
    
    try:
        with open(file_full_path, 'r') as infile:
            for line in infile:
                splitted = line.rstrip("\n").split("\t")
                path_info = splitted[-1].split("|")
                first_last = (path_info[0], path_info[-1])
                couple_shortest_path_dict[first_last].update(path_info)
    except FileNotFoundError:
        print(f"Error: File '{file_full_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return couple_shortest_path_dict


def track_shortest_paths(file_path, file_name, gene_pairs):
    """
    Tracks the list of shortest paths for each gene couple.

    Parameters:
    file_path (str): The directory path where the file is located.
    file_name (str): The name of the PathLinker output file.
    gene_pairs (set): A set of gene pairs.

    Returns:
    dict: A dictionary where keys are gene pairs and values are lists of shortest paths.
    """
    couple_list_of_shortest_paths_dict = {couple: [] for couple in gene_pairs}
    file_full_path = os.path.join(file_path, file_name)
    
    try:
        with open(file_full_path, 'r') as infile:
            for line in infile:
                splitted = line.rstrip("\n").split("\t")
                path_info = splitted[-1].split("|")
                first_last = (path_info[0], path_info[-1])
                couple_list_of_shortest_paths_dict[first_last].append(path_info)
    except FileNotFoundError:
        print(f"Error: File '{file_full_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return couple_list_of_shortest_paths_dict


def perform_enrichment_analysis(couple_shortest_path_dict):
    """
    Performs enrichment analysis for each gene couple.

    Parameters:
    couple_shortest_path_dict (dict): A dictionary where keys are gene pairs and values are sets of shortest path genes.
    """
    names = gseapy.get_library_name()
    count = 0

    for couple in couple_shortest_path_dict.keys():
        if couple != ('path', 'path'):
            gene1, gene2 = couple
            gene_list = list(couple_shortest_path_dict[couple])
            try:
                gseapy.enrichr(gene_list=gene_list, gene_sets='KEGG_2019_Human', outdir=f"Data/ENRICHR_GeneSetEnrichmentAnalysis_{gene1}_{gene2}")
            except ValueError:
                print(couple)
                count += 1
            except KeyError:
                print(couple)
                count += 1
    print(f"These couples give error, number: {count}")


def extract_signaling_pathways(couple_shortest_path_dict):
    """
    Extracts signaling pathways for each gene couple and writes the information to a file.

    Parameters:
    couple_shortest_path_dict (dict): A dictionary where keys are gene pairs and values are sets of shortest path genes.
    """
    output_file = os.path.join('Data', "EnrichRAllDiffGenePairs_CouplesShortestPathsSignalingPathway_Info.txt")
    with open(output_file, "a") as outfile:
        outfile.write("Gene1\tGene2\tSizeshortesPath\tDataset\tSignalingPathway\tMatchInPathway\tAdjP\tMinusLog10P\tMatchingGenesInPathway\n")

    for couple in couple_shortest_path_dict.keys():
        gene1, gene2 = couple
        try:
            report_file = os.path.join('Data', f"ENRICHR_GeneSetEnrichmentAnalysis_{gene1}_{gene2}/KEGG_2019_Human.human.enrichr.reports.txt")
            with open(report_file, 'r') as infile:
                infile.readline()  # Skip the header line
                for line in infile:
                    splitted = line.rstrip("\n").split("\t")
                    dataset = splitted[0]
                    pathway = splitted[1]
                    match = splitted[2]
                    adj_p = splitted[4]
                    log_p = -(np.log10(float(adj_p)))
                    if "signaling" in pathway:
                        with open(output_file, "a") as outfile:
                            outfile.write(f"{gene1}\t{gene2}\t{len(couple_shortest_path_dict[couple])}\t{dataset}\t{pathway}\t{match}\t{adj_p}\t{log_p}\t{splitted[-1]}\n")
        except FileNotFoundError:
            print(f"Error: Report file for couple {gene1}_{gene2} not found.")
        except Exception as e:
            print(f"An error occurred while processing {gene1}_{gene2}: {e}")


# if __name__ == "__main__":
#     data_folder = 'Data'
#     file_name = "paths.txt"
#     gene_pairs = get_gene_couples(data_folder, file_name)
#     couple_shortest_path_dict = lump_shortest_path_genes(data_folder, file_name, gene_pairs)
#     couple_list_of_shortest_paths_dict = track_shortest_paths(data_folder, file_name, gene_pairs)
#     perform_enrichment_analysis(couple_shortest_path_dict)
#     extract_signaling_pathways(couple_shortest_path_dict)

# data_folder = 'Data'
# file_name = "paths.txt"
# gene_pairs = get_gene_couples(data_folder, file_name)
# couple_shortest_path_dict = lump_shortest_path_genes(data_folder, file_name, gene_pairs)
# couple_list_of_shortest_paths_dict = track_shortest_paths(data_folder, file_name, gene_pairs)
# extract_signaling_pathways(couple_shortest_path_dict)
