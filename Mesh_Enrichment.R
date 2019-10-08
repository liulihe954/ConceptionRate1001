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

MESH_Enrichment_all3 = Parse_Results(Mesh_results_b_raw)
#raw_pvalue_distribution

#MESH_Enrichment_all3 = dplyr::filter(MESH_Enrichment_all3,pvalue <)


