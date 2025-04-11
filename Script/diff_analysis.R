diff_analysis = function(RNA.counts.normalized, feature, test, ref, pval = 0.05, file.name){
  
  ####TFs differential active 
  # Generating test and ref data
  test_i <- which(feature == test)
  ref_i <- which(feature == ref)
  
  mat_test <- as.matrix(RNA.counts.normalized[,test_i, drop = F])
  mat_ref <- as.matrix(RNA.counts.normalized[,ref_i, drop = F])
  
  # Generating NULL model (test, reference)
  dnull <- ttestNull(mat_test, mat_ref, per=1000, verbose = T)
  
  # Generating signature: Observed DEG scores
  signature <- rowTtest(mat_test, mat_ref) #Computes t-statistic for each gene
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic)) #Transform pvalues to z scores to standarize strength and direction
  signature <- na.omit(signature) #Remove NAs
  signature <- signature[,1] #Take scores as a numeric vector
  
  # Running msVIPER
  mra <- msviper(signature, net_regulons, dnull, verbose = F, minsize = 5)
  saveRDS(mra$es, file=paste0("Rerun/Monoculture/msViper_result_", file.name, ".rds"))
  
  
  # Plot DiffActive TFs
  pdf(paste0("Rerun/Monoculture/Differential_TFs_", file.name, ".pdf"))
  print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  dev.off()
  
  tfs_da = mra$es$p.value
  tfs_da = tfs_da[tfs_da < pval]
  tfs_names = names(tfs_da)
  
  nes = mra$es$nes
  nes = nes[tfs_names] #only significant
  up = names(nes[nes > 0]) #up active TFs
  down = names(nes[nes < 0]) #down active TFs
  
  #Map TFs targets
  targets_up = net$target[which(net$source%in%up)]
  targets_down = net$target[which(net$source%in%down)]
  common = intersect(targets_up, targets_down)
  if(length(common)!=0){
    targets_up = targets_up[-which(targets_up %in% common)]
    targets_down = targets_down[-which(targets_down %in% common)]    
  }
  
  ## Enrichment
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(RNA.counts.normalized), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  #Enrichment using targets of up regulated TFs
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = targets_up, columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  reac <- enrichPathway(gene    = entrz$ENTREZID,
                        organism     = 'human',
                        universe = universe$ENTREZID,
                        pvalueCutoff = 0.05)

  kegg <- enrichKEGG(gene    = entrz$ENTREZID,
                     organism     = 'human',
                     universe = universe$ENTREZID,
                     pvalueCutoff = 0.05)

  if(nrow(data.frame(kegg))!=0){
    pdf(paste0("Rerun/Monoculture/ORA_KEGG_UP_TFs_", file.name, ".pdf"))
    write.csv(kegg, paste0("Rerun/Monoculture/ORA_KEGG_UP_TFs_", file.name, ".csv"))
    print(dotplot(kegg))
    dev.off()
  }

  if(nrow(data.frame(reac))!=0){
    pdf(paste0("Rerun/Monoculture/ORA_Reactome_UP_TFs_", file.name, ".pdf"))
    write.csv(kegg, paste0("Rerun/Monoculture/ORA_REAC_UP_TFs_", file.name, ".csv"))
    print(dotplot(reac))
    dev.off()
  }

  #Enrichment using targets of down regulated TFs
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = targets_down, columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  reac <- enrichPathway(gene    = entrz$ENTREZID,
                        organism     = 'human',
                        universe = universe$ENTREZID,
                        pvalueCutoff = 0.05)

  kegg <- enrichKEGG(gene    = entrz$ENTREZID,
                     organism     = 'human',
                     universe = universe$ENTREZID,
                     pvalueCutoff = 0.05)

  if(nrow(data.frame(kegg))!=0){
    pdf(paste0("Rerun/Monoculture/ORA_KEGG_DOWN_TFs_", file.name, ".pdf"))
    write.csv(kegg, paste0("Rerun/Monoculture/ORA_KEGG_DOWN_TFs_", file.name, ".csv"))
    print(dotplot(kegg))
    dev.off()
  }

  if(nrow(data.frame(reac))!=0){
    pdf(paste0("Rerun/Monoculture/ORA_Reactome_DOWN_TFs_", file.name, ".pdf"))
    write.csv(kegg, paste0("Rerun/Monoculture/ORA_REAC_DOWN_TFs_", file.name, ".csv"))
    print(dotplot(reac))
    dev.off()
  }

  return(list(TFs_UP = up, TFs_DOWN = down, mra))
  
}
