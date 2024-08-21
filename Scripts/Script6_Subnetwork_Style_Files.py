#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 00:33:33 2024

@author: yavuzb2
"""
 
from Scripts import Script1_FunctionsForDataCuration as Datasets



def SubnetworkStyleFile(pair, PATH='Data/GenePairSpecificSubnetworks-CentralityMeasures/'):
    TFs = Datasets.TRRUST_Transcription_Factor_Set()
    OGs = Datasets.OncoKB_OG_TSG()['OG']
    TSGs = Datasets.OncoKB_OG_TSG()['TSG']
    RTKs = Datasets.RTKs()
    DrugTargets = Datasets.CMAP_drug_target_dictionary()
    sif_filename = PATH + 'GenePair_{}_{}_ShortestPathEnrichedPathways_Subnetwork_HIPPIE_ForCentralityMEasures.txt'.format(pair[0], pair[1])
    SubnetworkProteins = set()

    with open(sif_filename, "r") as infile:
        for line in infile:
            if ('Node' not in line) and ('nan' not in line):
                splitted = line.rstrip("\n").split("\t")
                SubnetworkProteins.add(splitted[0])
                SubnetworkProteins.add(splitted[1])
    
    # Prepare output file with TF, OG-TSG, RTK, etc.
    output_filename = PATH + 'GenePair_{}_{}_Subnetwork_Protein_Annotations.txt'.format(pair[0], pair[1])
    
    with open(output_filename, "w") as outfile:
        outfile.write("Protein\tTF\tOG-TSG\tRTK\tDrugTarget\n")
        
        for protein in SubnetworkProteins:
            is_tf = "TF" if protein in TFs else "notTF"
            
            if protein in OGs:
                og_tsg_status = "OG"
            elif protein in TSGs:
                og_tsg_status = "TSG"
            else:
                og_tsg_status = "neither"
                
            is_rtk = "RTK" if protein in RTKs else "notRTK"
            is_drug_target = "DrugTarget" if protein in DrugTargets.keys() else "notDrugTarget"
            
            outfile.write(f"{protein}\t{is_tf}\t{og_tsg_status}\t{is_rtk}\t{is_drug_target}\n")

    print(f"File saved: {output_filename}")

# Example usage:



if __name__ == "__main__":
    SubnetworkStyleFile(('ESR1', 'PIK3CA'))
    SubnetworkStyleFile(('BRAF', 'PIK3CA'))


########
###########3