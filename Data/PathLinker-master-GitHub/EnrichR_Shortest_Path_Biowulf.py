#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 07:54:42 2024

@author: yavuzb2
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 20:22:49 2024

@author: yavuzb2
"""

#get the list of gene couples with shortest path information from PathLinker
AllPairs = set()
with open("paths.txt","r") as infile:
    for line in infile:
        splitted = line.rstrip("\n").split("\t")
        PATH = splitted[-1].split("|")
        couple = (PATH[0],PATH[-1])
        AllPairs.add(couple)
        
###########################        
#Lump up all shortest path genes for each gene couple         
Couple_ShortestPathDict = {couple:set() for couple in AllPairs}

with open("paths.txt","r") as infile:
    for line in infile:
        splitted = line.rstrip("\n").split("\t")
        PATH = splitted[-1].split("|")
        FirstLast = (PATH[0],PATH[-1])
        Couple_ShortestPathDict[FirstLast] = Couple_ShortestPathDict[FirstLast].union(set(splitted[-1].split("|")))
                
###########################        
#Also keep track of list of shortest paths for each couple     
Couple_ListOfShortesPathsDict = {couple:[] for couple in AllPairs}

with open("paths.txt","r") as infile:
    for line in infile:
        splitted = line.rstrip("\n").split("\t")
        PATH = splitted[-1].split("|")
        FirstLast = (PATH[0],PATH[-1])
        Couple_ListOfShortesPathsDict[FirstLast].append(PATH)
                
##########
##################        
        
        
        
import gseapy        

names = gseapy.get_library_name()
#######3
count=0
for couple in Couple_ShortestPathDict.keys():
   
        gene1 = couple[0]
        gene2 = couple[1]
        Gene_List = list(Couple_ShortestPathDict[couple])
        try:
            gseapy.enrichr(gene_list=Gene_List, gene_sets='KEGG_2019_Human',outdir = "ENRICHR_GeneSetEnrichmentAnalysis_{}_{}".format(gene1,gene2))
        except ValueError:
            print(couple)
            count+=1

print("These couples gives error_ number: {}".format(count))
#!pip install gseapy
ShortestPathProteins = Couple_ShortestPathDict
# if you have conda (MacOS_x86-64 and Linux only)
#$ conda install -c bioconda gseapy
# Windows and MacOS_ARM64(M1/2-Chip)
#$ pip install gseapy
# assign dataframe, and use enrichr library data set 'KEGG_2016'
# assign a list object to enrichr
# gl = ['SCARA3', 'LOC100044683', 'CMBL', 'CLIC6', 'IL13RA1', 'TACSTD2', 'DKKL1', 'CSF1',
#      'SYNPO2L', 'TINAGL1', 'PTX3', 'BGN', 'HERC1', 'EFNA1', 'CIB2', 'PMP22', 'TMEM173']
with open("EnrichRAllDiffGenePairs_CouplesShortestPathsSignalingPathway_Info.txt","a") as outfile:
    outfile.write("Gene1"+"\t"+ "Gene2"+"\t"+ "SizeshortesPath"+"\t"+ "Dataset"+"\t"+ "SignalingPathway"+"\t"+ "MatchInPathway"+"\t"+ "AdjP"+"\t"+ "MinusLog10P"+"\t"+'MatchingGenesInPAthway'+"\n")

import numpy as np
for couple in Couple_ShortestPathDict.keys():
    gene1 =couple[0]
    gene2 = couple[1]
    with open("ENRICHR_GeneSetEnrichmentAnalysis_{}_{}/KEGG_2019_Human.human.enrichr.reports.txt".format(gene1,gene2)) as infile:
        first_line = infile.readline()
        for line in infile:
            #print(line)
            splitted = line.rstrip("\n").split("\t")
            dataset = splitted[0]
            pathway = splitted[1]
            match = splitted[2]
            adj_p = splitted[4]
            log_p = -(np.log10(float(adj_p )))
            if "signaling" in pathway:
                with open("EnrichRAllDiffGenePairs_CouplesShortestPathsSignalingPathway_Info.txt","a") as outfile:
                    outfile.write(gene1+"\t"+ gene2+"\t"+ str(len(ShortestPathProteins[couple]))+"\t"+ dataset+"\t"+ pathway+"\t"+ str(match)+"\t"+ str(adj_p) +"\t"+ str(log_p)+"\t"+ splitted[-1]+"\n")
                
            
 
            
            
            
