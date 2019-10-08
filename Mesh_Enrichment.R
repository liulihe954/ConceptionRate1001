##Function calls
source("Interpro_Enrichment_Founctions.R") 
## Data pre

# Read in total genesstr(Total_gene)
Total_gene_all = read.table("total.genes.slope.all.txt",sep = "") 
Total_gene_all$x = as.character(Total_gene_all$x)

# Sig genes
Sig_gene_all = read.table("sig.genes.slope.all.txt",sep = "") 
Sig_gene_all$Var1 = as.character(Sig_gene_all$Var1)

#
# convert ensemble to Entrez
# keys return the keys for the database contained in the MeSHdb object

key.symbol = keys(org.Bt.eg.db,  keytype = c("ENSEMBL"))
entrezUniverse = AnnotationDbi::select(org.Bt.eg.db, as.character(key.symbol), 
                        columns = c("ENTREZID"),keytype = "ENSEMBL") %>% 
  dplyr::distinct(ENSEMBL,.keep_all= TRUE)
dim(entrezUniverse)
# 
tmp1 = dplyr::left_join(Total_gene_all,entrezUniverse, by = c('x' = "ENSEMBL"))
Sig_list_out_entrez = list(all3_sig = unique(na.omit(tmp1$ENTREZID)))
tmp2 = dplyr::left_join(Sig_gene_all,entrezUniverse, by = c('Var1' = "ENSEMBL"))
Total_list_out_entrez = list(all3_total = unique(na.omit(tmp2$ENTREZID)))

#str(Sig_list_out_entrez)
#str(Total_list_out_entrez)

Total_gene_list_all3 = Sig_list_out_entrez
Sig_gene_list_all3 = Total_list_out_entrez
TestingSubsetNames = "All_three_slope_combined"

test = MESH_Enrich(Total_gene_list_all3,
                   Sig_gene_list_all3 ,
                   TestingSubsetNames,
                   Meshthres = 0.05,
                   MeshCate = c("D","G"),
                   #biomart="ensembl",
                   dataset="MeSH.Bta.eg.db",
                   #dataset= "btaurus_gene_ensembl",
                   #Identifier = "external_gene_name",
                   #attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                   keyword = "MESH_Enrichment_all3_1007")

load("MESH_Enrichment_all3_1007.RData")
#MESH_Enrichment_all3 = Parse_Results(Mesh_results_b_raw)


#MESH_Enrichment_all3 = dplyr::filter(MESH_Enrichment_all3,pvalue <)
Total_gene_list_all3 = Sig_list_out_entrez
Sig_gene_list_all3 = Total_list_out_entrez
TestingSubsetNames = "All_three_slope_combined"

# input pre
Total_gene_list_all3 = unique(unlist(Total_gene_list_all3));attributes(Total_gene_list_all3) = NULL
Sig_gene_list_all3 = unique(unlist(Sig_gene_list_all3));attributes(Sig_gene_list_all3) = NULL

### use funtion
meshParams <- new("MeSHHyperGParams", 
                  geneIds = Sig_gene_list_all3, 
                  universeGeneIds = Total_gene_list_all3, 
                  annotation = "MeSH.Bta.eg.db",
                  category = c("D"), database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
#str(meshParams)
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH_Term_ID", "MeSH_Term_Name", 
                  "Total_Genes", "DE_Genes", "P-value")
print(unique(out), row.names = F)

