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


####
GO_Enrich_Results_thres005_1014 = Go_Enrich_Plot(Total_gene_list,
                                                 Sig_gene_list,
                                                 TestingSubsetNames,GOthres = 0.05,
                                                 keyword = "GO_Enrichment_allslope_1015")


# formating
biomart="ensembl";dataset="btaurus_gene_ensembl";attributes = c("go_id","namespace_1003")
database = useMart(biomart);genome = useDataset(dataset, mart = database);gene = getBM(attributes,mart = genome)
namespace_index = dplyr::filter(gene,go_id != "",namespace_1003 != "")

#
load("GO_Enrichment_allslope_1015.RData")

#Interpro_Enrich_Results_thres005_1 = Parse_Interpro_Results(Interpro_results_b[[1]])
# slope 2
Interpro_Enrich_Results_thres005_1 = Parse_Results(GO_results_b[1], keyword= "slope1")
#names(Interpro_Enrich_Results_thres005_1)
Interpro_Enrich_Results_thres005_1 = 
  dplyr::select(Interpro_Enrich_Results_thres005_1,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue) %>% dplyr::left_join( namespace_index, by =c("GOID" ="go_id")) %>% dplyr::rename(go_id = GOID)
#
Interpro_Enrich_Results_thres005_2 = Parse_Results(GO_results_b[2], keyword= "slope2")
Interpro_Enrich_Results_thres005_2 = 
  dplyr::select(Interpro_Enrich_Results_thres005_2,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue)%>% dplyr::left_join( namespace_index, by =c("GOID" ="go_id")) %>% dplyr::rename(go_id = GOID)

# slope 3
Interpro_Enrich_Results_thres005_3 = Parse_Results(GO_results_b[3],keyword= "slope3")
Interpro_Enrich_Results_thres005_3 = 
  dplyr::select(Interpro_Enrich_Results_thres005_3,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue) %>% dplyr::left_join( namespace_index, by =c("GOID" ="go_id")) %>% dplyr::rename(go_id = GOID)

Interpro_Enrich_Results_thres005_4 = Parse_Results(GO_results_b[4], keyword= "slopeall")
Interpro_Enrich_Results_thres005_4 = 
  dplyr::select(Interpro_Enrich_Results_thres005_4,-ExternalLoss_total,-InternalLoss_sig) %>% 
  dplyr::arrange(pvalue) %>% dplyr::left_join( namespace_index, by =c("GOID" ="go_id")) %>% dplyr::rename(go_id = GOID)

# Output
require(openxlsx)
Interpro_Results <- list("Slope1" = Interpro_Enrich_Results_thres005_1, 
                         "Slope2" = Interpro_Enrich_Results_thres005_2,
                         "Slope3" = Interpro_Enrich_Results_thres005_3,
                         "Slope4" = Interpro_Enrich_Results_thres005_4)
write.xlsx(Interpro_Results,file = "GO_Results_all_1015.xlsx")

#######
ReduceDim_GO_Plot(Interpro_Enrich_Results_thres005_1,GOthres = 0.05, 
                  label_sizeCC = 0.4,
                  label_sizeBP = 0.4,
                  label_sizeMF = 0.4,Dataset_Name = "Slop1_sematic_sim_GO")
ReduceDim_GO_Plot(Interpro_Enrich_Results_thres005_2,GOthres = 0.05,
                  label_sizeCC = 0.4,
                  label_sizeBP = 0.4,
                  label_sizeMF = 0.4,Dataset_Name = "Slop2_sematic_sim_GO")
ReduceDim_GO_Plot(Interpro_Enrich_Results_thres005_3,GOthres = 0.05,
                  label_sizeCC = 0.4,
                  label_sizeBP = 0.4,
                  label_sizeMF = 0.4,Dataset_Name = "Slop3_sematic_sim_GO")

ReduceDim_GO_Plot(Interpro_Enrich_Results_thres005_4,GOthres = 0.05,
                  label_sizeCC = 0.4,
                  label_sizeBP = 0.3,
                  label_sizeMF = 0.2,Dataset_Name = "Slop4_sematic_sim_GOh")

#"GO:0038023""GO:0004970","GO:0005509","GO:0038023"
