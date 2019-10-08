suppressMessages(library(org.Bt.eg.db))
options(warn=-1); suppressMessages(library(meshr)); options(warn=0)
suppressMessages(library(MeSH.db))
suppressMessages(library(MeSH.Bta.eg.db))

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
#dim(entrezUniverse)
#
tmp1 = dplyr::left_join(Total_gene_all,entrezUniverse, by = c('x' = "ENSEMBL"))
Sig_list_out_entrez = list(all3_sig = unique(na.omit(tmp1$ENTREZID)))
tmp2 = dplyr::left_join(Sig_gene_all,entrezUniverse, by = c('Var1' = "ENSEMBL"))
Total_list_out_entrez = list(all3_total = unique(na.omit(tmp2$ENTREZID)))



##### second way
## FULL GENES

genes.back = data.frame(Total_gene_all$x)
colnames(genes.back) <- "ENSEMBL"
str(genes.back)
geneID.back <- merge(genes.back, entrezUniverse, by ="ENSEMBL")
str(geneID.back)
geneID2.back <- geneID.back[ !duplicated(geneID.back[,2]),]
str(geneID2.back)

## SIGNIFICANT GENES

genes.sig = data.frame(Sig_gene_all$Var1)
colnames(genes.sig) <- "ENSEMBL"
str(genes.sig)
geneID.sig <- merge(genes.sig, entrezUniverse, by ="ENSEMBL")
str(geneID.sig)
geneID2.sig <- geneID.sig[ !duplicated(geneID.sig[,2]),]
str(geneID2.sig)

##

ns = length(geneID2.sig[,1])
nt = length(geneID2.back[,1])
cat(paste("Significant Genes:", ns, " and Backgroung Genes:", nt - ns), "\n")

#######################

### Chemicals and Drugs
meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,2],
                  universeGeneIds = geneID2.back[,2],
                  annotation = "MeSH.Bta.eg.db",
                  category = "D",
                  database = "gene2pubmed",
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM,
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name",
                  "NT.Genes", "DE.Genes", "P-value")
print(unique(out), row.names = F)

### Phenomena and Processes

#table(duplicated(na.omit(tmp1$ENTREZID)))

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


str(Total_gene_list_all3)
str(Sig_gene_list_all3)

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


##
meshParams_g <- new("MeSHHyperGParams", 
                  geneIds = Sig_gene_list_all3, 
                  universeGeneIds = Total_gene_list_all3, 
                  annotation = "MeSH.Bta.eg.db",
                  category = c("G"), database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")


meshR_g <- meshHyperGTest(meshParams_g)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH_Term_ID", "MeSH_Term_Name", 
                  "Total_Genes", "DE_Genes", "P-value")
print(unique(out), row.names = F)
