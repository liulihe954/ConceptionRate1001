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


#str(Total_gene_list)
#str(Sig_gene_list)

  ##################################
  ##         Converting          ##
  #################################
# convert ensemble to Entrez
# keys return the keys for the database contained in the MeSHdb object
key.symbol = AnnotationDbi::keys(org.Bt.eg.db,  keytype = c("ENSEMBL"))
entrezUniverse = AnnotationDbi::select(org.Bt.eg.db, as.character(key.symbol), 
                                       columns = c("ENTREZID"),keytype = "ENSEMBL") %>% 
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)
#dim(entrezUniverse)
sig_genes_all = Sig_gene_list
total_genes_all = Total_gene_list

#library(limma)
Sig_list_out = list()
Total_list_out = list()

for (i in c(1:4)){
  ## for sig
  #i = 3
  tmp1 = data.frame(ENSEMBL = unlist(sig_genes_all[[i]]))
  #head(tmp1)
  tmp1 = dplyr::left_join(tmp1 ,entrezUniverse, by = c("ENSEMBL" = "ENSEMBL"))
  Sig_list_out[[i]] = tmp1;names(Sig_list_out)[i] = names(sig_genes_all)[i]
  ## for total 
  tmp2 = data.frame(ENSEMBL = unlist(total_genes_all[[i]]))
  tmp2 = dplyr::left_join(tmp2,entrezUniverse,by = c("ENSEMBL" = "ENSEMBL"))
  Total_list_out[[i]] = tmp2;names(Total_list_out)[i] = names(total_genes_all)[i]
}
##
# Keep only the entrez ID: then we have one vector for each element of the list (some format as always)
Sig_list_out_entrez = list()
Total_list_out_entrez = list()
for (i in c(1:4)){
  Sig_list_out_entrez[[i]] = unique(na.omit(data.frame(Sig_list_out[[i]])$ENTREZID))
  names(Sig_list_out_entrez)[i] = names(Sig_list_out)[i]
  Total_list_out_entrez[[i]] = unique(na.omit(data.frame(Total_list_out[[i]])$ENTREZID))
  names(Total_list_out_entrez)[i] = names(Total_list_out)[i]
}

#str(Sig_list_out_entrez)
#str(Total_list_out_entrez)
## save the convert for next time, it will take some time if run again
save(Total_list_out,Sig_list_out,Sig_list_out_entrez,Total_list_out_entrez,file = "ConvertName2Entrez.RData")
load("ConvertName2Entrez.RData")

#########################################################################################################


#########################################################################################################
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_lowest_path_bt,10)
# all_path
NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = 
  dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,
                ReactomeID = V2,
                Reactome_Description = V4, 
                Source = V5, 
                Species = V6)
#head(NCBI2Reactome_all_path_bt)
# all_react
NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = 
  dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V4,V6,V2,V3,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V4, 
                Reaction_Description = V6,
                ProteinID = V2,
                Protein_Description = V3,
                Source = V7, Species = V8)

#str(Sig_list_out)
#head(NCBI2Reactome_all_react_bt,50)
# turn data input as charactor
NCBI2Reactome_all_react_bt[] <-   lapply(NCBI2Reactome_all_react_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_lowest_path_bt[] <- lapply(NCBI2Reactome_lowest_path_bt, function(x) if(is.factor(x)) as.character(x) else x)
NCBI2Reactome_all_path_bt[] <-   lapply(NCBI2Reactome_all_path_bt, function(x) if(is.factor(x)) as.character(x) else x)

#########################################################################################################################
sig_genes_all = Sig_list_out_entrez
total_genes_all = Total_list_out_entrez
#InputSource = NCBI2Reactome_all_path_bt
TestingSubsetNames = "all 3 slopes combined"
#
Reactome_Enrichment_all3slope_1008 = 
  Reactome_Enrich(total_genes_all,
                  sig_genes_all,
                  Sig_list_out = Sig_list_out,
                  TestingSubsetNames,
                  NCBI2Reactome_all_path_bt,
                  Reacthres = 0.05,
                  keyword = "Reactome_Enrichment_all3slope_1009_all_path")
#
Reactome_Enrichment_all3slope_1008 = 
  Reactome_Enrich(total_genes_all,
                  sig_genes_all,
                  Sig_list_out = Sig_list_out,
                  TestingSubsetNames,
                  NCBI2Reactome_lowest_path_bt,
                  Reacthres = 0.05,
                  keyword = "Reactome_Enrichment_all3slope_1009_lowest_path")
#
Reactome_Enrichment_all3slope_1008 = 
  Reactome_Enrich(total_genes_all,
                  sig_genes_all,
                  Sig_list_out = Sig_list_out,
                  TestingSubsetNames,
                  NCBI2Reactome_all_react_bt,
                  Reacthres = 0.05,
                  keyword = "Reactome_Enrichment_all3slope_1009_all_react")


#load("Reactome_Enrichment_all3slope_1008.RData")

## Data pre


#Total_gene_list_sep
#Sig_gene_list_sep
# Looping index (if multiple, then multi-loop)
#TestingSubsetNames = names(Sig_gene_list)

