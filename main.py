#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 21:52:33 2024

@author: yavuzb2
"""

# main.py

import sys
#Replace PATH with the path of the folder CombinationTherapyTargets
sys.path.append('PATH/CombinationTherapyTargets')

# main.py

import CombinationTherapyTargets
print(CombinationTherapyTargets.__file__)


def run_data_curation_and_analysis_for_combination_therapy_targets1():
    import os

    import pandas as pd
    # Example usage of the functions
    from Scripts import Script1_FunctionsForDataCuration as S1
    gene_pairs = S1.ListOfGenePairsWithCoexisitngMutations()
    S1.HIPPIE_PPI_Network(url='https://www.ndexbio.org/viewer/networks/89dd3925-3718-11e9-9f06-0ac135e8bacf')
    kegg_pathway_gene_dict = S1.KEGG2019_Pathway_Gene_Dictionary()
    transcription_factors = S1.TRRUST_Transcription_Factor_Set()
    oncogenes_and_tumor_suppressors = S1.OncoKB_OG_TSG()
    rtks = S1.RTKs()
    
    # Print results or perform further processing
    print(f"Gene pairs: {gene_pairs}")
    print(f"KEGG Pathways: {kegg_pathway_gene_dict}")
    print(f"Transcription Factors: {transcription_factors}")
    print(f"Oncogenes and Tumor Suppressors: {oncogenes_and_tumor_suppressors}")
    print(f"RTKs: {rtks}")

#run_data_curation_and_analysis_for_combination_therapy_targets1()



def run_data_curation_and_analysis_for_combination_therapy_targets3():
    import os

    import pandas as pd
    from Scripts import Script3_ShortestPathGenesPathwayEnrichment_EnrichR as S3
    data_folder = 'Data'
    file_name = "paths.txt"
    gene_pairs = S3.get_gene_couples(data_folder, file_name)
    couple_shortest_path_dict = S3.lump_shortest_path_genes(data_folder, file_name, gene_pairs)
    couple_list_of_shortest_paths_dict = S3.track_shortest_paths(data_folder, file_name, gene_pairs)
    S3.perform_enrichment_analysis(couple_shortest_path_dict)
    S3.extract_signaling_pathways(couple_shortest_path_dict)

#run_data_curation_and_analysis_for_combination_therapy_targets3()

def run_data_curation_and_analysis_for_combination_therapy_targets4():
    import os

    import pandas as pd
    from Scripts import Script4_Processing_EnrichR_Results_EnrichedPathways as S4
    gene_pair_pathways = S4.process_enrichr_results()
    S4.write_signaling_pathways_to_file()
    

#run_data_curation_and_analysis_for_combination_therapy_targets4()



   


def run_data_curation_and_analysis_for_combination_therapy_targets5():
    import os

    import pandas as pd
    # # Read gene pair pathways Script 5, S5
    from Scripts import Script5_GenePair_PageRankShortestPaths as S5
    GenePairDict_OnePathways = S5.read_gene_pair_pathways()[0]
    GenePairDict_TwoPathways = S5.read_gene_pair_pathways()[1]
    GenePairDict_MorethanTwoPathways = S5.read_gene_pair_pathways()[2]
    gene_couple_pathway_shortest_path_proteins = S5.gene_couple_shortest_path_genes()
    S5.create_folder_in_data(folder_name = 'GenePairSpecificSubnetworks-CentralityMeasures')

    Merged_GenePairDict_OneTwo_MoreThanTwoPathways = {}
    Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_OnePathways)
    Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_TwoPathways)
    Merged_GenePairDict_OneTwo_MoreThanTwoPathways.update(GenePairDict_MorethanTwoPathways)
    
    
    for gene_pair in Merged_GenePairDict_OneTwo_MoreThanTwoPathways.keys():
        S5.gene_pair_specific_subnetwork(gene_pair,gene_pair_dict_pathway_info = Merged_GenePairDict_OneTwo_MoreThanTwoPathways, network_file =  "Data/HIPPI_PPI_network_self_loops_removed.txt")   
    

    for pair in Merged_GenePairDict_OneTwo_MoreThanTwoPathways.keys():
        Targets = S5.SelectTargetsWithBetweennessCentrality(pair,PATH = 'Data/GenePairSpecificSubnetworks-CentralityMeasures/')
        target = '+'.join(list(Targets))
        with open('Data/GenePair_DrugTargetsFromBetweennessCentrality_IQR_GreaterThan_Q3.txt','a') as outfile:
            outfile.write(pair[0]+"\t"+pair[1]+"\t"+target+'\n')
    # Example usage
    filename = 'Data/GenePair_DrugTargetsFromBetweennessCentrality_IQR_GreaterThan_Q3.txt'
    intersection_nodes = S5.process_gene_pair_interactions(GenePairDict_OnePathways, 
                                                        GenePairDict_TwoPathways, 
                                                        GenePairDict_MorethanTwoPathways, 
                                                      filename)
    print(intersection_nodes)

#run_data_curation_and_analysis_for_combination_therapy_targets5()

def run_data_curation_and_analysis_for_combination_therapy_targets6():
    import os

    import pandas as pd
    from Scripts import Script6_Subnetwork_Style_Files as S6
    #please add your gene pairs
    GenePairsForSubnetwork = [('ESR1', 'PIK3CA'), ('BRAF', 'PIK3CA')]
    for pair in GenePairsForSubnetwork:
        S6.SubnetworkStyleFile(pair)

#run_data_curation_and_analysis_for_combination_therapy_targets6()




if __name__ == "__main__":
    # Example usage of the functions
   run_data_curation_and_analysis_for_combination_therapy_targets1()
   run_data_curation_and_analysis_for_combination_therapy_targets3()
   run_data_curation_and_analysis_for_combination_therapy_targets4()
   run_data_curation_and_analysis_for_combination_therapy_targets5()
   run_data_curation_and_analysis_for_combination_therapy_targets6()
   
   
   
   
   
   
   
   
