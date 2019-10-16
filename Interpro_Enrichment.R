##Function calls
source("Interpro_Enrichment_Founctions.R") 
## Data pre

##################################
##    the combined slope       ##
#################################
# Read in total genesstr(Total_gene)
Total_gene_all = read.table("total.genes.slope.all.txt",sep = "") 
Total_gene_all = as.character(Total_gene_all$x)

# Sig genes
Sig_gene_all = read.table("sig.genes.slope.all.txt",sep = "") 
Sig_gene_all$Var1 = as.character(Sig_gene_all$Var1)
Sig_gene_all = as.character(Sig_gene_all$Var1)

##################################
##    the first three           ##
#################################
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
                       Slope3_T = Total_gene3$x,
                       Slope_all_T = Total_gene_all)

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
Sig_gene_list= list(Slope1_S = Sig_gene$x,
                    Slope2_S = Sig_gene2$x,
                    Slope3_S = Sig_gene3$x,
                    Slope_all_S = Sig_gene_all)

# Looping index (if multiple, then multi-loop)
TestingSubsetNames = names(Sig_gene_list)

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
                  keyword = "Interpro_Enrichment_allslope_1014")



load("Interpro_Enrichment_allslope_1014.RData")
## Massage the output
# match family annotation
library(tidyverse)
match_family = read.table("entry.list",sep = "\t",header = T)
match_family = as_tibble(match_family) %>% dplyr::select(ENTRY_AC,ENTRY_TYPE) %>% 
  rename(InterproID = ENTRY_AC,Type =ENTRY_TYPE)


# formating
load("Interpro_Enrichment_allslope_1014.RData")
#Interpro_Enrich_Results_thres005_1 = Parse_Interpro_Results(Interpro_results_b[[1]])
# slope 2
Interpro_Enrich_Results_thres005_1 = Parse_Results(Interpro_results_b[1], keyword= "slope1")
#names(Interpro_Enrich_Results_thres005_1)
Interpro_Enrich_Results_thres005_1 = 
  dplyr::select(Interpro_Enrich_Results_thres005_1,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue_r) %>% 
  dplyr::left_join(match_family,by = c("InterproID"="InterproID"))

Interpro_Enrich_Results_thres005_2 = Parse_Results(Interpro_results_b[2], keyword= "slope2")
Interpro_Enrich_Results_thres005_2 = 
  dplyr::select(Interpro_Enrich_Results_thres005_2,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue_r) %>% 
  dplyr::left_join(match_family,by = c("InterproID"="InterproID"))

# slope 3
Interpro_Enrich_Results_thres005_3 = Parse_Results(Interpro_results_b[3],keyword= "slope3")
Interpro_Enrich_Results_thres005_3 = 
  dplyr::select(Interpro_Enrich_Results_thres005_3,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue_r) %>% 
  dplyr::left_join(match_family,by = c("InterproID"="InterproID"))

Interpro_Enrich_Results_thres005_4 = Parse_Results(Interpro_results_b[4], keyword= "slopeall")
Interpro_Enrich_Results_thres005_4 = 
  dplyr::select(Interpro_Enrich_Results_thres005_4,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue_r) %>% 
  dplyr::left_join(match_family,by = c("InterproID"="InterproID"))

# Output
require(openxlsx)
Interpro_Results <- list("Slope1" = Interpro_Enrich_Results_thres005_1, 
                         "Slope2" = Interpro_Enrich_Results_thres005_2,
                         "Slope3" = Interpro_Enrich_Results_thres005_3,
                         "Slope4" = Interpro_Enrich_Results_thres005_4)
write.xlsx(Interpro_Results,file = "Interpro_Results_all_1015.xlsx")

