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
                  label_sizeBP = 0.4,
                  label_sizeMF = 0.4,Dataset_Name = "Slop4_sematic_sim_GOh")

#"GO:0038023""GO:0004970","GO:0005509","GO:0038023"

go_test = c("GO:0008066")
test =  AnnotationDbi::select(org.Bt.eg.db, keys = go_test, keytype="GO", "ENTREZID")

names(Interpro_Enrich_Results_thres005_1)
test_a = dplyr::filter(Interpro_Enrich_Results_thres005_2, go_id == "GO:0008066")
test_a 
ReduceDim_GO_Plot = function(Enrich_Out,
                             GOthres = 0.001,
                             label_sizeCC = 0.4,
                             label_sizeBP = 0.4,
                             label_sizeMF = 0.4,
                             Database = "org.Bt.eg.db",
                             measure="Jiang",combine=NULL,
                             Dataset_Name){
  # load libraries + download ref database
  library(GOSemSim);library(corrplot);library(tidyverse)
  do.call(library,list(Database))
  semData_BP <- godata(paste(Database), ont="BP", computeIC=T)
  semData_MF <- godata(paste(Database), ont="MF", computeIC=T)
  semData_CC <- godata(paste(Database), ont="CC", computeIC=T)
  # selection + formating: for each category we have one vector containing all the sig GO terms
  BP_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "biological_process") %>% 
    dplyr::select(go_id) %>% unlist();attributes(BP_List) = NULL # name is an attribute and we dont them, so set null
  CC_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "cellular_component") %>% 
    dplyr::select(go_id) %>% unlist();attributes(CC_List) = NULL
  MF_List = dplyr::filter(Enrich_Out,pvalue<=GOthres & namespace_1003 == "molecular_function") %>% 
    dplyr::select(go_id) %>% unlist();attributes(MF_List) = NULL
  ### Now we are trying to get all similarity matrix ready. N x N, symetric, diag = 1
  # For BP
  
  goSimMatrix_BP = GOSemSim::mgoSim(BP_List,
                                    BP_List,
                                    semData=semData_BP,measure=measure,combine = combine)
  suspectID_BP = rownames(goSimMatrix_BP)[is.na(goSimMatrix_BP[,1])]
  if (length(suspectID_BP) != 0){BP_List_new = setdiff(BP_List,suspectID_BP)
  message(length(suspectID_BP)," invalid ID captured in BP: ",suspectID_BP,", thus been removed!")
  } else {BP_List_new = BP_List;message("Nice! All IDs are valid in BP!")}
  goSimMatrix_BP_new = GOSemSim::mgoSim(BP_List_new,
                                        BP_List_new,
                                        semData=semData_BP,measure=measure,combine = combine)
  colnames(goSimMatrix_BP_new) = paste(BP_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)])
  rownames(goSimMatrix_BP_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% BP_List_new)],BP_List_new)
  # For CC
  goSimMatrix_CC = GOSemSim::mgoSim(CC_List,
                                    CC_List,
                                    semData=semData_CC,measure=measure,combine = combine)
  suspectID_CC = rownames(goSimMatrix_CC)[is.na(goSimMatrix_CC[,1])]
  if (length(suspectID_CC) != 0){CC_List_new = setdiff(CC_List,suspectID_CC)
  message(length(suspectID_CC)," invalid ID captured in CC: ",suspectID_CC,", thus been removed!")
  } else {CC_List_new = CC_List;message("Nice! All IDs are valid in CC!")}
  goSimMatrix_CC_new = GOSemSim::mgoSim(CC_List_new,
                                        CC_List_new,
                                        semData=semData_CC,measure=measure,combine =combine)
  colnames(goSimMatrix_CC_new) = paste(CC_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)])
  rownames(goSimMatrix_CC_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% CC_List_new)],CC_List_new)
  # For MF
  goSimMatrix_MF = GOSemSim::mgoSim(MF_List,
                                    MF_List,
                                    semData=semData_MF,measure=measure,combine = combine)
  suspectID_MF = rownames(goSimMatrix_MF)[is.na(goSimMatrix_MF[,1])]
  if (length(suspectID_MF) != 0){MF_List_new = setdiff(MF_List,suspectID_MF)
  message(length(suspectID_MF)," invalid ID captured in MF: ",suspectID_MF,", thus been removed!")
  } else {MF_List_new = MF_List;message("Nice! All IDs are valid in MF!")}
  goSimMatrix_MF_new = GOSemSim::mgoSim(MF_List_new,
                                        MF_List_new,
                                        semData=semData_MF,measure=measure,combine = combine)
  colnames(goSimMatrix_MF_new) = paste(MF_List_new,Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)])
  rownames(goSimMatrix_MF_new) = paste(Enrich_Out$GO_Name[(Enrich_Out$go_id %in% MF_List_new)],MF_List_new)
  # Now we take the results and plot
  pdf(paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".pdf",sep = ""))
  corrplot(goSimMatrix_CC_new,title = "Semantic_Similarity_Measure_CC",
           tl.col = "black", tl.cex = label_sizeCC, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_BP_new,title = "Semantic_Similarity_Measure_BP",
           tl.col = "black", tl.cex = label_sizeBP, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  corrplot(goSimMatrix_MF_new,title = "Semantic_Similarity_Measure_MF",
           tl.col = "black", tl.cex = label_sizeMF, 
           method = "shade", order = "hclust", 
           hclust.method = "centroid", is.corr = FALSE,mar=c(0,0,1,0))
  dev.off()
  message(dim(goSimMatrix_CC_new)[1],",",
          dim(goSimMatrix_BP_new)[1],",",
          dim(goSimMatrix_MF_new)[1]," GOs ploted in CC, BP and MF, respectively",
          ", cutting at, ",GOthres)
  require(openxlsx)
  CorMatrix <- list("CorMat_BP" = data.frame(goSimMatrix_BP), 
                    "CorMat_CC" = data.frame(goSimMatrix_CC),
                    "CorMat_MF" = data.frame(goSimMatrix_MF))
  write.xlsx(CorMatrix, row.names=TRUE,
             file = paste("Semantic_Similarity_Measure_",
                          Dataset_Name,"_",
                          formatC(GOthres, format = "e", digits = 0),".xlsx",sep = ""))
  save(goSimMatrix_CC,goSimMatrix_BP,goSimMatrix_MF,
       file = paste("Semantic_Similarity_Measure_",Dataset_Name,"_",formatC(GOthres, format = "e", digits = 0),".RData",sep = ""))
  message("Nice! Excels, Plots exported and RData saved!")
}