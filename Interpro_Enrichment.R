##Function calls
source("Interpro_Enrichment_Founctions.R") # Function Call
## Data pre
# Read in total genesstr(Total_gene)
Total_gene = read.table("total.genes.slope1.txt",sep = "") 
Total_gene$x = as.character(Total_gene$x)
#
Total_gene2 = read.table("total.genes.slope2.txt",sep = "") 
Total_gene2$x = as.character(Total_gene2$x)
#
Total_gene3 = read.table("total.genes.slope3.txt",sep = "") 
Total_gene3$x = as.character(Total_gene3$x)

Total_gene_list = list(Slope1_T = Total_gene$x,
                       Slope2_T = Total_gene2$x,
                       Slope3_T = Total_gene3$x)

# # Read in sig genes
Sig_gene = read.csv("sig.genes.slope1.txt",sep = "") 
Sig_gene$x = as.character(Sig_gene$x)
Sig_gene_list = list(Sig = Sig_gene$x)
#
Sig_gene2 = read.csv("sig.genes.slope2.txt",sep = "") 
Sig_gene2$x = as.character(Sig_gene2$x)
#
Sig_gene3 = read.csv("sig.genes.slope3.txt",sep = "") 
Sig_gene3$x = as.character(Sig_gene3$x)

Sig_gene_list = list(Slope1_S = Sig_gene$x,
                     Slope2_S = Sig_gene2$x,
                     Slope3_S = Sig_gene3$x)


# Looping index (if multiple, then multi-loop)
TestingSubsetNames = names(Sig_gene_list)
TestingSubsetNames = TestingSubsetNames[-1]
Sig_gene_list = Sig_gene_list[2:3]
Total_gene_list = Total_gene_list[2:3]

############################
#@    Enrichment loops   ##
##########################
Interpro_Enrich_Results_thres005 = 
  InterPro_Enrich(Total_gene_list,
                  Sig_gene_list,
                  TestingSubsetNames,
                  IPthres = 0.05,
                  biomart="ensembl",
                  dataset="btaurus_gene_ensembl",
                  Identifier = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                  keyword = "Interpro_Enrichment_other2slopes_1007")

load("Interpro_Enrichment_other2slopes_1007.RData")
## Massage the output
# match family annotation
library(tidyverse)
match_family = read.table("entry.list",sep = "\t",header = T)
match_family = as_tibble(match_family) %>% dplyr::select(ENTRY_AC,ENTRY_TYPE) %>% 
  rename(InterproID = ENTRY_AC,Type =ENTRY_TYPE)


# formating
load("Interpro_Enrichment_other2slopes_1007.RData")
#Interpro_Enrich_Results_thres005_1 = Parse_Interpro_Results(Interpro_results_b[[1]])
# slope 2
Interpro_Enrich_Results_thres005_2 = Parse_Results(Interpro_results_b[1], keyword= "slope2")
Interpro_Enrich_Results_thres005_2 = 
  dplyr::select(Interpro_Enrich_Results_thres005_2,ID,Description,Total_gene,Significant_gene,pvalue, HitPerc ) %>% 
  dplyr::arrange(pvalue) %>% 
  dplyr::left_join(match_family,by = c("ID"="InterproID"))

# slope 3
Interpro_Enrich_Results_thres005_3 = Parse_Results(Interpro_results_b[2],keyword= "slope3")
Interpro_Enrich_Results_thres005_3 = 
  dplyr::select(Interpro_Enrich_Results_thres005_3,ID,Description,Total_gene,Significant_gene,pvalue, HitPerc ) %>% 
  dplyr::arrange(pvalue) %>% 
  dplyr::left_join(match_family,by = c("ID"="InterproID"))

# Output
write.csv(Interpro_Enrich_Results_thres005_2,"Interpro_Enrich_Results_thres005_2.csv",row.names = F) 
write.csv(Interpro_Enrich_Results_thres005_3,"Interpro_Enrich_Results_thres005_3.csv",row.names = F) 

#write.csv(Interpro_Enrich_Results_thres005_final,"Interpro_Enrich_Results_thres005_final.csv",row.names = F) 

##===============================#
## second round - combined     ##
##=============================#
## Data pre
# Read in total genesstr(Total_gene)
Total_gene_all = read.table("unique_all_slopes_Total_genes.txt",sep = "") 
Total_gene_all$x = as.character(Total_gene_all$x)
# Sig genes
Sig_gene_all = read.table("genes_repeat_2_3times.txt",sep = "") 
Sig_gene_all$Var1 = as.character(Sig_gene_all$Var1)
#
Total_gene_list_all3 = list(Total_gene_all = Total_gene_all$x)
Sig_gene_list_all3 = list(Sig_gene_all = Sig_gene_all$Var1)
TestingSubsetNames = "All_three_slope_combined"


# run all interpro
Interpro_Enrich_Results_thres005_all = 
  InterPro_Enrich(total_genes_all = Total_gene_list_all3,
                  sig_genes_all = Sig_gene_list_all3,
                  TestingSubsetNames,
                  IPthres = 0.05,
                  biomart="ensembl",
                  dataset="btaurus_gene_ensembl",
                  Identifier = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                  keyword = "Interpro_Enrichment_all_three_combined")

# match info
library(tidyverse)
match_family = read.table("entry.list",sep = "\t",header = T)
match_family = as_tibble(match_family) %>% dplyr::select(ENTRY_AC,ENTRY_TYPE) %>% 
  rename(InterproID = ENTRY_AC,Type =ENTRY_TYPE)


# formating
load("Interpro_Enrichment_all_three_combined.RData")
#Interpro_Enrich_Results_thres005_1 = Parse_Interpro_Results(Interpro_results_b[[1]])
# slope 2
Interpro_Enrich_Results_thres005_all = Parse_Results(Interpro_results_b, keyword= "all_three_slope_combined")
Interpro_Enrich_Results_thres005_all = 
  dplyr::select(Interpro_Enrich_Results_thres005_all,ID,Description,Total_gene,Significant_gene,pvalue, HitPerc) %>% 
  dplyr::arrange(pvalue) %>% 
  dplyr::left_join(match_family,by = c("ID"="InterproID"))

# Output
write.csv(Interpro_Enrich_Results_thres005_all,"Interpro_Enrich_Results_thres005_all.csv",row.names = F) 
