library(biomaRt)
library(org.Bt.eg.db)

# pick a target
go_test = c("GO:0008066")







library(biomaRt)
library(org.Bt.eg.db)
# get info from biomart
database = useMart("ensembl")
genome = useDataset("btaurus_gene_ensembl", mart = database)
gene = getBM(attributes = c("ensembl_gene_id","go_id","name_1006"),mart = genome)
# all the go from biomart
all_go1 = unique(na.omit(gene$go_id))[-1]
length(all_go1) # total 15118
# all the go from org.Bt.eg.db
all_go2 = AnnotationDbi::keys(org.Bt.eg.db,keytype = c("GO"))
length(all_go2) # total 9032
# intersect
table(all_go2 %in% all_go1)
table(all_go1 %in% all_go2)

# test one of them
head(all_go1[!(all_go1 %in% all_go2)])

sessionInfo()










# check the target in biomart
test1 = dplyr::filter(gene,go_id == go_test)
# Check the target in org.Bt.eg.db
test2 =  AnnotationDbi::select(org.Bt.eg.db, keys = go_test, keytype="GO", "ENSEMBL")


# all the go from biomart
all_go1 = unique(na.omit(gene$go_id))[-1]
length(all_go1)
# all the go from org.Bt.eg.db
all_go2 = AnnotationDbi::keys(org.Bt.eg.db,keytype = c("GO"))
length(all_go2)


# intersect
table(all_go2 %in% all_go1)
table(all_go1 %in% all_go2)

# biomart gives more
# common - only 8716