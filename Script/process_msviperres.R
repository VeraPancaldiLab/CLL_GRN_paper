
library(pheatmap)
library(igraph)

#--------------------------------------------
# Load and classify msVIPER result files
#--------------------------------------------

auto_dir <- "Results/msVIPER/Autologous/Objects/"
mono_dir <- "Results/msVIPER/Monoculture/Objects/"

# Separate by experimental condition
# Extract only files that contain msVIPER results
autofiles <- list.files(path = auto_dir, pattern = '^msViper.*\\.rds$', full.names = T)
monofiles <- list.files(path = mono_dir, pattern = '^msViper.*\\.rds$', full.names = T)

files <- c(autofiles, monofiles)

files = autofiles

#--------------------------------------------
# Parse msVIPER results, filter significant TFs
#--------------------------------------------
d <- list()     # Store filtered results
tfs <- list()   # Store significant TF names
nes <- list()   # Store NES values for each TF

for (f in files) {
  ddf <- as.data.frame(readRDS(f)[-2])                        # Read RDS and drop 2nd item
  d[[f]] <- ddf[ddf$p.value < 0.05, ]                         # Filter by p-value
  tfs[[f]] <- rownames(d[[f]])                                # Save TF names
  nes[[f]] <- d[[f]]$nes                                      # Save NES values
  names(nes[[f]]) <- tfs[[f]]                                 # Name NES by TFs
}

#--------------------------------------------
# Create unified list of unique TFs
#--------------------------------------------
tfsuniqueauto <- unique(unlist(tfs[autofiles]))
tfsuniquemono <- unique(unlist(tfs[monofiles]))
tfsunique <- unique(unlist(tfs[files]))

#--------------------------------------------
# Build NES matrix across all samples
#--------------------------------------------
mat <- matrix(0, nrow=length(tfsunique), ncol=length(files))
rownames(mat) <- tfsunique
colnames(mat) <- files

for (f in files) {
  mat[tfs[[f]], f] <- nes[[f]]
}

# Save NES matrix
write.table(mat, 'TFannot.txt', quote=FALSE, sep='\t')

#--------------------------------------------
# Prepare ordered file list for AUTO condition
#--------------------------------------------
filesordauto <- autofiles[c(2,3,4,1, 6,7,8,5,10,11,12,9)]

# Build AUTO NES matrix
mata <- matrix(0, nrow=length(tfsuniqueauto), ncol=length(filesordauto))
rownames(mata) <- tfsuniqueauto
colnames(mata) <- filesordauto

for (f in filesordauto) {
  mata[tfs[[f]], f] <- nes[[f]]
}

# Split AUTO by patient
matap1 <- mata[, grep('Patient1', colnames(mata))]
matap1 <- matap1[rowSums(matap1) > 0, ]
matap2 <- mata[, grep('Patient2', colnames(mata))]
matap2 <- matap2[rowSums(matap2) > 0, ]
matap3 <- mata[, grep('Patient3', colnames(mata))]
matap3 <- matap3[rowSums(matap3) > 0, ]

#--------------------------------------------
# Prepare ordered file list for MONO condition
#--------------------------------------------
filesordmono <- monofiles[c(2,3,4,1, 6,7,8,5,10,11,12,9)]

# Build MONO NES matrix
matm <- matrix(0, nrow=length(tfsuniquemono), ncol=length(filesordmono))
rownames(matm) <- tfsuniquemono
colnames(matm) <- filesordmono

for (f in filesordmono) {
  matm[tfs[[f]], f] <- nes[[f]]
}

# Split MONO by patient
matmp1 <- matm[, grep('Patient1', colnames(matm))]
matmp1 <- matmp1[rowSums(matmp1) > 0, ]
matmp2 <- matm[, grep('Patient2', colnames(matm))]
matmp2 <- matmp2[rowSums(matmp2) > 0, ]
matmp3 <- matm[, grep('Patient3', colnames(matm))]
matmp3 <- matmp3[rowSums(matmp3) > 0, ]

#--------------------------------------------
# Generate heatmaps for visualization
#--------------------------------------------
pdf('tfheat-clus.pdf', width=10, height=15)

# MONO heatmaps
pheatmap(matm, cexRow=0.5, main='All mono', scale='none')
pheatmap(matmp1, cexRow=0.8, cluster_cols=FALSE, main='Patient1 mono', scale='none')
pheatmap(matmp2, cexRow=0.8, cluster_cols=FALSE, main='Patient2 mono', scale='none')
pheatmap(matmp3, cexRow=0.8, cluster_cols=FALSE, main='Patient3 mono', scale='none')

# AUTO heatmaps
pheatmap(mata, cexRow=0.5, main='All auto', scale='none')
pheatmap(matap1, cexRow=0.8, cluster_cols=FALSE, main='Patient1 auto', scale='none')
pheatmap(matap2, cexRow=0.8, cluster_cols=FALSE, main='Patient2 auto', scale='none')
pheatmap(matap3, cexRow=0.8, cluster_cols=FALSE, main='Patient3 auto', scale='none')

dev.off()

# Full dataset heatmap (AUTO + MONO)
pdf('tfheat-clusall.pdf', width=10, height=15)
pheatmap(mat, cexRow=0.5, main='All', scale='none', cluster_cols=FALSE)
dev.off()

#--------------------------------------------
# Combine AUTO and MONO matrices
#--------------------------------------------
matall <- rbind(mata, matm)

# Correlation between samples (conditions)
allcondcor <- cor(t(matall))
corcond <- cor(mat)
cortf <- cor(t(mat))

# Create TF co-correlation network

Gcortf <- graph_from_adjacency_matrix(cortf, weight=TRUE, mode='undirected')
write.table(
  cbind(get.edgelist(Gcortf), E(Gcortf)$weight),
  'Gcortfnet.txt',
  quote=FALSE, row.names=FALSE
)

#--------------------------------------------
# Annotate positive/negative TFs per patient
#--------------------------------------------
p1tfsannot <- apply(matp1, 2, function(x) which(x > 0))
p2tfsannot <- apply(matp2, 2, function(x) which(x > 0))
p3tfsannot <- apply(matp3, 2, function(x) which(x > 0))

p1tfsnegannot <- apply(matp1, 2, function(x) which(x < 0))
p2tfsnegannot <- apply(matp2, 2, function(x) which(x < 0))
p3tfsnegannot <- apply(matp3, 2, function(x) which(x < 0))
