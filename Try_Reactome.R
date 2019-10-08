## Data pre
# Read in total genesstr(Total_gene)
Total_gene_all = read.table("total.genes.slope.all.txt",sep = "") 
Total_gene_all$x = as.character(Total_gene_all$x)
# Sig genes
Sig_gene_all = read.table("sig.genes.slope.all.txt",sep = "") 
Sig_gene_all$Var1 = as.character(Sig_gene_all$Var1)
str(Sig_gene_all)

# convert ensemble to Entrez
# keys return the keys for the database contained in the MeSHdb object
key.symbol = keys(org.Bt.eg.db,  keytype = c("ENSEMBL"))
entrezUniverse = AnnotationDbi::select(org.Bt.eg.db, as.character(key.symbol), 
                                       columns = c("ENTREZID"),keytype = "ENSEMBL") %>% 
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)
#dim(entrezUniverse)
#
tmp1 = dplyr::left_join(Sig_gene_all,entrezUniverse, by = c('Var1' = "ENSEMBL"))
Sig_list_out_entrez = list(all3_sig = unique(na.omit(tmp1$ENTREZID)))
tmp2 = dplyr::left_join(Total_gene_all,entrezUniverse, by = c('x' = "ENSEMBL"))
Total_list_out_entrez = list(all3_total = unique(na.omit(tmp2$ENTREZID)))

# Read in database
NCBI2Reactome_lowest_path = read.csv("NCBI2Reactome.txt",sep = "\t",header = F)
NCBI2Reactome_lowest_path_bt = dplyr::filter(NCBI2Reactome_lowest_path, V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  dplyr::rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_lowest_path_bt,10)

NCBI2Reactome_all_path = read.csv("NCBI2Reactome_All_Levels.txt",sep = "\t",header = F)
NCBI2Reactome_all_path_bt = dplyr::filter(NCBI2Reactome_all_path,V6 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V4,V5,V6) %>% 
  rename(EntrezID = V1,ReactomeID = V2,Reactome_Description = V4, Source = V5,Species = V6)
#head(NCBI2Reactome_all_path_bt)

NCBI2Reactome_all_react = read.csv("NCBI2Reactome_PE_Reactions.txt",sep = "\t",header = F)
NCBI2Reactome_all_react_bt = dplyr::filter(NCBI2Reactome_all_react,V8 == "Bos taurus") %>% 
  dplyr::select(V1,V2,V3,V4,V6,V7,V8) %>% 
  dplyr::rename(EntrezID = V1,ProteinID = V2,Protein_Description = V3,ReactionID = V4, 
                Reaction_Description = V6,Source = V7, Species = V8)
#head(NCBI2Reactome_all_react_bt)

#########################################################################################################################
sig_genes_all = Sig_list_out_entrez
total_genes_all = Total_list_out_entrez
#InputSource = NCBI2Reactome_all_path_bt
TestingSubsetNames = "all 3 slopes combined"
#
Reactome_Enrichment_all3slope_1008 = 
  Reactome_Enrich(total_genes_all,
                  sig_genes_all,
                  TestingSubsetNames,
                  NCBI2Reactome_all_path_bt,
                  Reacthres = 0.05,
                  keyword = "Reactome_Enrichment_all3slope_1008_all_path")
#
Reactome_Enrichment_all3slope_1008 = 
  Reactome_Enrich(total_genes_all,
                  sig_genes_all,
                  TestingSubsetNames,
                  NCBI2Reactome_lowest_path_bt,
                  Reacthres = 0.05,
                  keyword = "Reactome_Enrichment_all3slope_1008_lowest_path")
#
Reactome_Enrichment_all3slope_1008 = 
  Reactome_Enrich(total_genes_all,
                  sig_genes_all,
                  TestingSubsetNames,
                  NCBI2Reactome_all_react_bt,
                  Reacthres = 0.05,
                  keyword = "Reactome_Enrichment_all3slope_1008_all_react")


#load("Reactome_Enrichment_all3slope_1008.RData")


