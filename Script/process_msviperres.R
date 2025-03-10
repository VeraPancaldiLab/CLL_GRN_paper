files=list.files()

files=files[grep('msViper', files)]
autofiles=files[grep('^msViper', files)]


monofiles=files[grep('monomsViper', files)]

d=list()
tfs=list()
nes=list()
for (f in files[grep('msViper',files)]){
ddf=as.data.frame(readRDS(f)[-2])

d[[f]]=ddf[which(ddf[,'p.value']<0.05),]
tfs[[f]]=rownames(d[[f]])
nes[[f]]=d[[f]]$nes
names(nes[[f]])=tfs[[f]]
}



d=list()
tfs=list()
nes=list()
for (f in files[grep]){
ddf=as.data.frame(readRDS(f)[-2])

d[[f]]=ddf[which(ddf[,'p.value']<0.05),]
tfs[[f]]=rownames(d[[f]])
nes[[f]]=d[[f]]$nes
names(nes[[f]])=tfs[[f]]
}





##make a list of unique TFs

tfsuniqueauto=unique(unlist(tfs[autofiles]))


tfsuniquemono=unique(unlist(tfs[monofiles]))

tfsunique=unique(unlist(tfs[files]))

mat=matrix(0, nrow=length(tfsunique), ncol=length(files))
colnames(mat)=files
rownames(mat)=tfsunique

for (f in files){
mat[tfs[[f]],f]=nes[[f]]
}
allcondcor=cor(t(matall))


 ##calculate correlation between conditions
corcond=cor(mat)

cortf=cor(t(mat))


 Gcortf=graph_from_adjacency_matrix(cortf, weight=T, mode='undirected')
write.table(cbind(get.edgelist(Gcortf), E(Gcortf)$weight), 'Gcortfnet.txt', quote=F, row.names=F)


write.table(mat, 'TFannot.txt', quote=F, sep='\t')






####consider auto
##list conditions in the right time order and separate patients
filesordauto=autofiles[c(2,3,4,1, 6,7,8,5,10,11,12,9)]

##make full heatmap
mata=matrix(0, nrow=length(tfsuniqueauto), ncol=length(filesordauto))
colnames(mata)=filesordauto
rownames(mata)=tfsuniqueauto

for (f in filesordauto){
mata[tfs[[f]],f]=nes[[f]]
}

matap1=mata[,grep('Patient1', colnames(mata))]
matap1=matap1[which(rowSums(matap1)>0),]
matap2=mata[,grep('Patient2', colnames(mata))]
matap2=matap2[which(rowSums(matap2)>0),]
matap3=mata[,grep('Patient3', colnames(mata))]
matap3=matap3[which(rowSums(matap3)>0),]
 


####consider mono
##list conditions in the right time order and separate patients
filesordmono=monofiles[c(2,3,4,1, 6,7,8,5,10,11,12,9)]

##make full heatmap
matm=matrix(0, nrow=length(tfsuniquemono), ncol=length(filesordmono))
colnames(matm)=filesordmono
rownames(matm)=tfsuniquemono

for (f in filesordmono){
matm[tfs[[f]],f]=nes[[f]]
}

matmp1=matm[,grep('Patient1', colnames(matm))]
matmp1=matmp1[which(rowSums(matmp1)>0),]
matmp2=matm[,grep('Patient2', colnames(matm))]
matmp2=matmp2[which(rowSums(matmp2)>0),]
matmp3=matm[,grep('Patient3', colnames(matm))]
matmp3=matmp3[which(rowSums(matmp3)>0),]


 pdf('tfheat-clus.pdf', width=10, height=15)
 pheatmap(matm, cexRow=0.5, main='All mono', scale='none')
 pheatmap(matmp1, cexRow=0.8, cluster_cols=F, main='Patient1 mono', scale='none')
 pheatmap(matmp2, cexRow=0.8,cluster_cols=F, main='Patient2 mono', scale='none')
 pheatmap(matmp3, cexRow=0.8,cluster_cols=F, main='Patient3 mono', scale='none')

 pheatmap(mata, cexRow=0.5,  main='All auto', scale='none')
 pheatmap(matap1, cexRow=0.8, cluster_cols=F, main='Patient1 auto', scale='none')
 pheatmap(matap2, cexRow=0.8,cluster_cols=F, main='Patient2 auto', scale='none')
 pheatmap(matap3, cexRow=0.8,cluster_cols=F, main='Patient3 auto', scale='none')

dev.off()

pdf('tfheat-clusall.pdf', width=10, height=15)
pheatmap(mat,cexRow=0.5, main='All', scale='none', cluster_cols=F)
dev.off()


matall=cbind(mata, matm)

allcondcor=cor(t(matall))


 ##calculate correlation between conditions
corcond=cor(mat)

cortf=cor(t(mat))


 Gcortf=graph_from_adjacency_matrix(cortf, weight=T, mode='undirected')
write.table(cbind(get.edgelist(Gcortf), E(Gcortf)$weight), 'Gcortfnet.txt', quote=F, row.names=F)

##find which TFs are associated with which point

p1tfsannot=apply(matp1,2, function(x){return(which(x>0))})
p2tfsannot=apply(matp2,2, function(x){return(which(x>0))})
p3tfsannot=apply(matp3,2, function(x){return(which(x>0))})

p1tfsnegannot=apply(matp1,2, function(x){return(which(x<0))})
p2tfsnegannot=apply(matp2,2, function(x){return(which(x<0))})
p3tfsnegannot=apply(matp3,2, function(x){return(which(x<0))})