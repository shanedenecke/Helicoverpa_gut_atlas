#!/usr/bin/env Rscript

####################################################  Import Libraries #################################################### 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(edgeR))
shhh(library(gplots))
shhh(library(Mfuzz))

rep.mean=function(data,col){
  sample=select(data,matches(col))
  return(rowMeans(sample))
}

column.clean=function(x){
    x1=gsub('d_','G_Artificial_',x)
    x2=gsub('p_','G_Plant_',x1)
    x3=gsub('(L[0-9])_prot','\\1_MG_prot',x2)
    return(x3)
}
#################################################### Set working directory  #################################################### 

setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/fuzzyClustering',showWarnings = F)

####################################################  Import data  #################################################### 
trans.raw=fread('./outputs/omics_clean/Ha_proteomics_transcriptomics_Replicates.csv')
trans.sum=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')


####################################################  Prepare Heatmap data #################################################### 
colnames(trans.raw)=column.clean(colnames(trans.raw))
trans.raw=select(trans.raw,ha_geneid,(matches('Plant',perl=T)))
trans.matrix=as.matrix(trans.raw,rownames = 'ha_geneid') 
trans.filter=trans.matrix[apply(trans.matrix,1,sd)>1,] 
y <- DGEList(counts = trans.filter, group=rep(1:as.integer(dim(trans.matrix)[2]/4),each=4))
y <- calcNormFactors(y)
z <- cpm(y, normalized.lib.size=TRUE) 
scaledata <- t(scale(t(z))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.

####################################################  Graph Heatmap data #################################################### 
tiff('./outputs/fuzzyClustering/Heatmap.tiff')
hm=heatmap.2(z,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none")
dev.off()

####################################################  Prepare Fuzzy data #################################################### 

set.seed(12346)
colnames(trans.sum)=column.clean(colnames(trans.sum))
fuzzy.matrix=trans.sum %>% select(ha_geneid,matches('Plant.+mean')) %>% 
  .[complete.cases(.)] %>% filter(!duplicated(ha_geneid)) %>% as.matrix(rownames = 'ha_geneid') 
fuzzy.filter=fuzzy.matrix[apply(fuzzy.matrix,1,sum)>1,] 

### Make expression dataset
gut.exp=ExpressionSet(assayData=fuzzy.filter)
gut.exp.f=filter.std(gut.exp,min.std=0,visu=F)
gut.exp.s <- standardise(gut.exp.f)
m1 <- mestimate(gut.exp.s)


####################################################  Make Dmins Plot #################################################### 
clust.step=Dmin(gut.exp.s, m=m1, crange=seq(2,10,1), repeats=3,visu = F) %>% 
  data.table(`Minimum Centroid Distance`=.,`Number of Clusters`=2:10)

gp=ggplot(data=clust.step,aes(y=`Minimum Centroid Distance`,x=`Number of Clusters`))
gp=gp+geom_bar(stat='identity')
gp=gp+scale_x_continuous(breaks=2:10)
gp=gp+theme_bw()
gp=gp+theme(axis.title=element_text(size=20,face='bold'),
            panel.grid=element_blank(),
            axis.text.x=element_text(size=14,face='bold'),
            axis.text.y=element_text(size=14,face='bold'),
            rect=element_rect(fill=NULL),
            legend.text=element_text(size=18),
            legend.title=element_text(size=26,face='bold'),
            strip.text = element_text(size=20,face='bold'),
            plot.title = element_text(size=28,face='bold',hjust = 0.5))
ggsave(plot=gp,filename='./outputs/fuzzyClustering/Custer_size_step.tiff',device='tiff',height=8,width=14)


####################################################  Make Fuzzy Plot #################################################### 


c.8 <- mfuzz(gut.exp.s,c=4,m=m1)
#c.8 <- mfuzz(gut.exp.s,c=6,m=m1)
#c.8 <- mfuzz(gut.exp.s,c=8,m=m1)
par(mar=c(1,1,1,1))
tiff('./outputs/fuzzyClustering/L5_fuzzy.tiff')
fuzzy.plot=mfuzz.plot2(gut.exp.s,cl=c.8,mfrow=c(2,2),
                       time.labels = c('FG','AMG','MMG','PMG','HG'),
                       xlab='Gut Compartment',x11=F,min.mem=.7,centre=T)
print(fuzzy.plot)
dev.off()


####################################################  Extract Fuzzy Genes #################################################### 
## annotate best cluster data and get best member
cluster=data.table(best.clust=c.8[[3]],ha_geneid=names(c.8[[3]]))
membership=c.8[[4]] %>% data.table()

for(i in 1:length(colnames(membership))){colnames(membership)[i]=paste0('A',colnames(membership)[i])} ### make column numeric
membership=membership %>% rowwise() %>% mutate(maximum=max(A1,A2,A3,A4)) %>% data.table() 
membership=cbind(cluster,membership)
fin.membership=membership %>% filter(maximum>.7) %>% arrange(desc(best.clust)) %>% data.table() 
mem.tab=fin.membership$best.clust %>% table()


####################################################  Annotate Fuzzy Clusters #################################################### 
key=data.table(Gut_specific=1263,AMG_up=680,HG_specific=2507,FG_specific=1623) %>% t() ### for .7 and 4 clusters

ml=list()
for(i in 1:4){
  n=dim(fin.membership[best.clust==i])[1]
  ml[[i]]=data.table(num=i,name=names(key[,][key[,]>n-10 & key[,]<n+10]))
}
final.fuzzy.key=rbindlist(ml)

cl_name=c()
for(i in fin.membership$best.clust){cl_name=c(cl_name,final.fuzzy.key[num==i]$name)}
fin.membership$cluster_name=cl_name
write.xlsx(fin.membership,'./outputs/Gene_Subset_expression/TableS8_Fuzzy_Clustering.xlsx',row.names = F)


### take a peak which pH genes come up
#ph=fread('./inputs/acid_base/HelArm_pH_genes.csv') %>% filter(Category!='SLC_2') %>% mutate(ha_geneid=paste0('LOC',ha_geneid)) %>%
#merge(fin.membership,by='ha_geneid') %>% arrange(cluster_name)

####################################################  Write Fuzzy Gene groups #################################################### 

for(i in unique(fin.membership$cluster_name)){
  fil=fin.membership %>% filter(cluster_name==i)
  writeLines(as.character(fil$ha_geneid),paste0('./outputs/gene_groups/Cluster_',as.character(i),'_Gene_List.txt'))
}




