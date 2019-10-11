##Function calls
source("Interpro_Enrichment_Founctions.R") 
suppressMessages(library(org.Bt.eg.db))
options(warn=-1); suppressMessages(library(meshr)); options(warn=0)
suppressMessages(library(MeSH.db))
suppressMessages(library(MeSH.Bta.eg.db))
## Data pre

# raw data for retrive MESHid and all details linked
#KEY = keys(MeSH.db, keytype = "MESHID")
#List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
#Match_List = dplyr::select(List, MESHID, MESHTERM)
# head(Match_List) 

# Prepare Bta database
#key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
#list_Bta = select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
#  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY == c("D","G")) %>% 
#  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
# head(list_Bta,30)
keyword_outer = "MeshDB"
DB = paste(keyword_outer,".RData",sep = "")
load(DB)

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
#####=============================================================#######
#####=============================================================#######
Total_gene_list_all3 = Total_list_out_entrez
Sig_gene_list_all3 = Sig_list_out_entrez
#str(Total_gene_list_all3)
#str(Sig_gene_list_all3)
TestingSubsetNames = names(total_genes_all)

Mesh_Enrichment_all3slope_1010 = 
  MESH_Enrich(Total_gene_list_all3,
              Sig_gene_list_all3,
              TestingSubsetNames,
              Sig_list_out = Sig_list_out,
              Meshthres = 0.05,
              MeshCate = c("D","G"),
              dataset="MeSH.Bta.eg.db",
              keyword = "Mesh_Enrichment_allslope_1011_G")
                   

#####=============================================================#######
#####=============================================================#######
##############################
### formating the results  ##
##############################
All_Results_List = c("Mesh_Enrichment_allslope_1011_G.RData")
All_Keywords_List = c("Mesh_Enrich_Slope_final_005_1011_G.xlsx")
# loop for outputs
for (i in seq_along(All_Results_List)){
  tmp_results_name = All_Results_List[i]
  tmp_key_name = All_Keywords_List[i]
  load(tmp_results_name)
  # set up selecting index
  compile_select_index = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","findG","hitsPerc")
 # slope1
  slope1 = Parse_Results(Mesh_results_b[1]) 
  names(slope1) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  slope1 = dplyr::select(slope1,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  # slope2
  slope2 = Parse_Results(Mesh_results_b[2]) 
  names(slope2) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  slope2 = dplyr::select(slope2,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  # slope3
  slope3 = Parse_Results(Mesh_results_b[3]) 
  names(slope3) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  slope3 = dplyr::select(slope3,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  # slope_all
  slopeall = Parse_Results(Mesh_results_b[4]) 
  names(slopeall) = c("MeshID","MeshTerm","Total_Genes","Significant_Genes","pvalue","ExternalLoss_total","InternalLoss_sig","findG","hitsPerc")
  slopeall = dplyr::select(slopeall,compile_select_index)  #%>% dplyr::left_join(match_family,by=c("InterproID" = "InterproID"))
  #
  Slope_final <- list("slope1" = slope1,
                      "slope2" = slope2,
                       "slope3" = slope3,
                       "slope4" = slopeall)
  require(openxlsx)
  write.xlsx(Slope_final,file = tmp_key_name)
}




