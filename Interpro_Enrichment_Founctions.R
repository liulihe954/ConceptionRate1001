library(org.Bt.eg.db)
library(meshr)
library(MeSH.db)
library(MeSH.Bta.eg.db)
library(tidyverse)
#########################################################################################################################
#########################################################################################################################
InterPro_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           IPthres = 0.05,
                           biomart="ensembl",
                           dataset="btaurus_gene_ensembl",
                           #host="http://www.ensembl.org",
                           Identifier = "external_gene_name",
                           attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Interpro_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Interpro_results_b = list()
  Interpro_results_b_raw = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr);library(tidyverse)# load pkg
  ## GetInterpro : bosTau 
  database = useMart(biomart)
  genome = useDataset(dataset, mart = database)
  gene = getBM(attributes,mart = genome)
  ##
  InterproName = unique(gene[,c("interpro","interpro_description")]) %>% arrange(interpro)
  Interpro = na.omit(InterproName$interpro)[-1]
  Name = na.omit(InterproName$interpro_description)[-1]
  #
  if (Identifier == "ensembl_gene_id"){genesInterpro = unique(subset(gene,interpro != "")$ensembl_gene_id)
  } else if (Identifier == "external_gene_name") {
    genesInterpro= unique(subset(gene,interpro != "")$external_gene_name);
    genesInterpro = genesInterpro[-1]
  } else {message("Sorry, we only have ensembel and names available as identifier, please use one of the followings: 
                  ensembl_gene_id OR external_gene_name.")}
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Interpro domains to check: ",length(Interpro)," with total number of names: ",length(Name))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesInterpro])
    S = length(sig.genes[sig.genes %in% genesInterpro]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(Interpro=character(),
                     Name=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(Interpro)){
      if (j%%1000 == 0) {message("tryingd on Interpro ",j," - ",Interpro[j]," - ",Name[j])}
      if (Identifier == "ensembl_gene_id"){
        gENEs = subset(gene, interpro == Interpro[j])$ensembl_gene_id
      } else if (Identifier == "external_gene_name") {
        gENEs = subset(gene, interpro == Interpro[j])$external_gene_name
      }
      m = length(total.genes[total.genes %in% gENEs]) # genes from target interpro and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target interpro also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(Interpro = Interpro[j], 
                       Name = Name[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Interpro_results_b_raw[[i]] = final_raw; names(Interpro_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Interpro raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= IPthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("InterproID","Interpro_Name", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Interpro_results_b[[i]] = final;names(Interpro_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched Interpro")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Interpro = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Interpro_results_b, Interpro_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Interpro domains found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",IPthres)
  message("Nice! - Interpro enrichment finished and data saved")}
#########################################################################################################################
Parse_Results = function(Results_List,keyword = "Which D.B"){
  all_enrich = data.frame(ID=character(),
                          Description=character(),
                          Total_gene=numeric(),
                          Significant_gene=numeric(),
                          pvalue=numeric(),
                          ExternalLoss_total = character(),
                          InternalLoss_total = character(),
                          HitPerc = numeric(),
                          stringsAsFactors=FALSE)
  for (i in 1:length(Results_List)){
    len = dim(data.frame(Results_List[i]))[1]
    if (len> 0){
      tmp = data.frame(Results_List[i])
      names(tmp) = names(all_enrich)
      all_enrich = rbind(all_enrich,tmp)
    }
  }
  #all_enrich_KEGG <- all_enrich_KEGG %>% dplyr::group_by(KEGG.ID) %>% dplyr::distinct()
  total_hits = dim(all_enrich)[1]
  total_modules = length(Results_List)
  print(paste("In database: ",keyword,"-",total_hits,"hits found in",total_modules,"tested modules: ",names(Results_List)))
  return(ParseResults = all_enrich)
}
#########################################################################################################################
#########################################################################################################################

#####=============================================================#######
MESH_Enrich = function(total_genes_all,
                       sig_genes_all,
                       TestingSubsetNames,
                       Meshthres = 0.05,
                       MeshCate = c("D","G"),
                       #biomart="ensembl",
                       dataset="MeSH.Bta.eg.db",
                       #dataset= "btaurus_gene_ensembl",
                       #Identifier = "external_gene_name",
                       #attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                       keyword = "MESH_Enrichment"){
  #total.genes = Total_list_out_entrez_test
  #sig.genes = Sig_list_out_entrez_test
  #TestingSubsetNames = TestingSubsetNames_test
  total_enrich = 0                        
  raw_pvalue_all = numeric()
  Mesh_results_b = list()
  Mesh_results_b_raw = list()
  library(MeSH.db);library(MeSH.Bta.eg.db);library(tidyverse);library(gage);library(magrittr)
  library(ggplot2);library(biomaRt);library(MeSH.Bta.eg.db)
  #========================================================================#
  ### raw data for retrive MESHid and all details linked
  #
  #KEY = keys(MeSH.db, keytype = "MESHID")
  #List = select(MeSH.db, keys = KEY, columns = columns(MeSH.db), keytype = "MESHID")
  ##List = select(MeSH.db, keys = KEY[1:3], columns = columns(MeSH.db), keytype = "MESHID")
  #Match_List = dplyr::select(List, MESHID, MESHTERM)
  ##head(Match_List) 
  ### Prepare Bta database
  #
  #key_Bta <- keys(MeSH.Bta.eg.db, keytype = "MESHID")
  #list_Bta = MeSHDbi::select(MeSH.Bta.eg.db, keys = key_Bta, columns = columns(MeSH.Bta.eg.db)[-4], keytype = "MESHID") %>% 
  #  dplyr::select(GENEID,MESHCATEGORY,MESHID,SOURCEID) %>% dplyr::filter(MESHCATEGORY == MeshCate) %>% 
  #  dplyr::left_join(Match_List,by= c("MESHID" = "MESHID"))
  #========================================================================#
  # alternatively
  keyword_outer = "MeshDB"
  DB = paste(keyword_outer,".RData",sep = "")
  load(DB)
  #Sig_list_out_entrez_test2
  #Total_list_out_entrez_test2
  # Get index
  genesMesh = unique(list_Bta$GENEID)
  MeshID = unique(list_Bta$MESHID)#;MeshID = MeshID[1:1000] # delete
  MeshTerm = unique(list_Bta$MESHTERM)
  #head(unique(MeshID),200)
  #length(genesGO)
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Mesh to check: ",length(MeshID)," with total number of names: ",length(MeshTerm))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% genesMesh])
    S = length(sig.genes[sig.genes %in% genesMesh]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(MeshID=character(),
                     MeshTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in c(1:length(MeshID))){
      if (j%%1000 == 0) {message("tryingd on MeshID ",j," - ",MeshID[j]," - ",MeshTerm[j])}
      gENEs = unique(subset(list_Bta, MESHID == MESHID[j])$GENEID)
      m = length(total.genes[total.genes %in% gENEs]) # genes from target  and in our dataset
      s = length(sig.genes[sig.genes %in% gENEs]) # # genes from target  also in the non-preserved module
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      #length(gENEs);Pval
      tmp = data.frame(MeshID= MeshID[j], 
                       MeshTerm = MeshTerm[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Mesh_results_b_raw[[i]] = final_raw; names(Mesh_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Mesh raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Meshthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("MeshID","MeshTerm", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Mesh_results_b[[i]] = final;names(Mesh_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched Mesh")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Mesh = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Mesh_results_b, Mesh_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant MeshIDs found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",Meshthres)
  message("Nice! - Mesh enrichment finished and data saved")}
#####=============================================================#######



#####=============================================================#######
Reactome_Enrich = function(total_genes_all,
                           sig_genes_all,
                           TestingSubsetNames,
                           InputSource,
                           Sig_list_out,
                           Reacthres = 0.05,
                           #biomart="ensembl",
                           #dataset="btaurus_gene_ensembl",
                           #host="http://www.ensembl.org",
                           #attributes = c("ensembl_gene_id","external_gene_name","interpro","interpro_description"),
                           keyword = "Reactome_Enrichment"){
  total_enrich = 0
  raw_pvalue_all = numeric()
  Reactome_results_b = list()
  Reactome_results_b_raw = list()
  library(ggplot2);library(biomaRt);library(gage);library(magrittr);library(tidyverse)# load pkg
  Reactome_gene =   unique(InputSource[,1])
  ReactomeID =      unique(InputSource[,2])
  ReactomeName =    unique(InputSource[,3])
  #ReactomeID = ReactomeID[1:300]
  message("Total Number of module/subsets to check: ",length(TestingSubsetNames))
  message("Total Number of Reactome to check: ",length(ReactomeID)," with total number of names: ",length(ReactomeName))
  #pdf(paste(trimws(keyword),".pdf",sep = ""))
  for (i in c(1:(length(TestingSubsetNames)))){
    #i = 1
    message("working on dataset #",i," - ",TestingSubsetNames[i])
    sig.genes = unlist(sig_genes_all[i]);attributes(sig.genes) = NULL
    total.genes = unlist(total_genes_all[i]);attributes(total.genes) = NULL
    #length(sig.genes)
    #length(total.genes)
    # total genes in the non-preserved module
    N = length(total.genes[total.genes %in% Reactome_gene])
    S = length(sig.genes[sig.genes %in% Reactome_gene]) #
    ExternalLoss_total = paste((length(total.genes) - N),round((length(total.genes) - N)/N,3),sep = "/")
    ExternalLoss_sig = paste((length(sig.genes) - S),round((length(sig.genes) - S)/S,3),sep = "/")
    out = data.frame(ReactomeID=character(),
                     ReactomeTerm=character(),
                     totalG=numeric(),
                     sigG=numeric(),
                     Pvalue=numeric(),
                     ExternalLoss_total = character(),
                     ExternalLoss_sig = character(),
                     findG =  character())
    message("Module size of ",TestingSubsetNames[i],": ", length(sig.genes))
    for(j in 1:length(ReactomeID)){
      # j = 101
      if (j%%100 == 0) {message("tryingd on Reactome ",j," - ",ReactomeID[j]," - ",ReactomeName[j])}
      target = ReactomeID[j]
      gENEs = unique(subset(InputSource, ReactomeID == target)$EntrezID)
      m = length(total.genes[total.genes %in% gENEs]) 
      findG = sig.genes[sig.genes %in% gENEs]
      s = length(findG)
      orig_list = data.frame(Sig_list_out[[i]]) %>% dplyr::filter(ENTREZID == findG)
      PastefindG = paste(orig_list$ENTREZID, collapse="/")
      M = matrix(c(s,S-s,m-s,N-m-S+s),byrow = 2, nrow = 2)
      Pval = round(fisher.test(M, alternative ="g")$p.value,100)
      tmp = data.frame(ReactomeID = ReactomeID[j], 
                       ReactomeName = ReactomeName[j], 
                       totalG = m, 
                       sigG = s, 
                       Pvalue = Pval, 
                       ExternalLoss_total = ExternalLoss_total,
                       ExternalLoss_sig = ExternalLoss_sig,
                       findG = PastefindG)
      out = rbind(out,tmp)}
    # put all palues in a box
    raw_pvalue_all = append(raw_pvalue_all,out$Pvalue,length(raw_pvalue_all))
    # raw complilation starts
    final_raw = out[order(out$Pvalue),];colnames(final_raw) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final_raw = final_raw %>% dplyr::mutate(hitsPerc = Significant_Genes*100 / Total_Genes)
    Reactome_results_b_raw[[i]] = final_raw; names(Reactome_results_b_raw)[i] = paste(TestingSubsetNames[i],"with",dim(final_raw)[1],"enriched Reactomeid raw")
    # raw complilation ends
    # selection starts - select those has 4 more gene in common and pvalue smaller than 0.05
    ot = subset(out,totalG > 4 & Pvalue <= Reacthres)
    final = ot[order(ot$Pvalue),];colnames(final) = c("ReactomeID","ReactomeName", "Total_Genes", "Significant_Genes", "pvalue_r","ExternalLoss_total","InternalLoss_sig","findG")
    final = final %>% mutate(hitsPerc = (Significant_Genes*100)/Total_Genes)
    Reactome_results_b[[i]] = final;names(Reactome_results_b)[i] = paste(TestingSubsetNames[i],"with",dim(final)[1],"enriched ReactomeID")
    # selection ends
    message("Significant Enrichment Hits:",nrow(final))
    total_enrich = total_enrich + nrow(final)
    ##
    #   print(final %>% 
    #           top_n(dim(final)[1], wt= -pvalue)%>% 
    #           ggplot(final, aes( x = hitsPerc,
    #                      y = GO_Name,
    #                      colour = pvalue,
    #                      size = Significant_Genes)) +
    #           geom_point() +
    #           theme_gray()+
    #          labs(title= paste("GO Enrichment in module",
    #                              TestingSubsetNames[i])), 
    #                              x="Hits (%)", y="GO term", 
    #                              colour="p value", size="Count")+
    #      theme(axis.text.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.text.y = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.x = element_text(size = 8,color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(axis.title.y = element_text(size = 8, color = "black",vjust = 0.5, hjust = 0.5))+
    #      theme(plot.title = element_text(size = 12,color = "black", face = "bold", vjust = 0.5, hjust = 0.5))
  }
  #  dev.off()
  raw_pvalue_index = seq(0.05,1,by=0.05)
  raw_pvalue_sum = numeric()
  for( z in seq_along(raw_pvalue_index)){raw_pvalue_sum[z] = length(which(raw_pvalue_all <= raw_pvalue_index[z]))}
  raw_pvalue_distribution = data.frame(index = raw_pvalue_index,counts_Reactome = raw_pvalue_sum)
  #raw_pvalue_distribution
  save(Reactome_results_b, Reactome_results_b_raw, raw_pvalue_distribution, file = paste(trimws(keyword),".RData",sep = ""))
  message(total_enrich," significant Reactome domains found within ",
          length(TestingSubsetNames)," modules/subsets", 
          " at the significance level of ",Reacthres)
  message("Nice! - Reactome enrichment finished and data saved")}

#####=============================================================#######

