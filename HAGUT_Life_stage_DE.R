#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(Mfuzz))
shhh(library(VennDiagram))
shhh(library(writexl))

######################################## import data ######################################## 
#setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
setwd('../../Dropbox/Shane/omics_projects/Ha_midgut_atlas/')
dir.create('./outputs/Gut_DEGs',showWarnings = F)
dir.create('./outputs/fuzzyClustering',showWarnings = F)

merged.omics=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')
key=fread('./inputs/keys/protein2gids.map',col.names = c('ha_prot_acc','ha_geneid','GeneID')) %>% select(ha_geneid,GeneID) %>% unique.data.frame()

###### Create Unigene protein key
#prot.key=key=fread('./inputs/keys/protein2gids.map',col.names = c('ha_prot_acc','ha_geneid','GeneID')) %>% select(ha_geneid,ha_prot_acc) %>% unique.data.frame()
a=fread('./inputs/keys/HelArm_NCBI_table.csv') %>% select(GeneID,`Protein product`,Length) 
colnames(a)=c('ha_geneid','ha_prot_acc','Length')
a$ha_geneid=paste0('LOC',a$ha_geneid) %>% gsub('LOC9','9',.)

ncbi.l=list()
for(i in unique(a$ha_geneid)){
  sub=a[ha_geneid==i]
  m=arrange(sub,desc(Length))[1]
  ncbi.l[[i]]=m
}
prot.key=rbindlist(ncbi.l)


trans.sum=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')


####################################### Define functions #################################
#x=fread('./inputs/DE_tables/L2MG__C2__L2-UP.csv')$sampleA
panos.rename=function(x){
  out=c()
  for(i in x){
    if(x %in% key$GeneID){
      out=c(out,paste0('LOC',key[GeneID %in% i]$ha_geneid))
    }else{
      print('test')
      out=c(out,NA)
    }
  }
  return(out)
}

#######################################  Import and clean DE tables ####################################### 
l2.gut=fread('./inputs/DE_tables/L2MG__C2__L2-UP.csv') %>% select(1,4,6,7) %>% rename(panos.id=sampleA,LogFC=logCPM,FDR=Geneid,PValue=FDR) %>% 
  mutate(ha_geneid=panos.rename(panos.id),Larval_Stage='L2') %>% filter(FDR<1e-4) %>% filter(ha_geneid!='LOC')
l3.gut=fread('./inputs/DE_tables/L3MG__C3__L3-UP.csv') %>% select(1,4,6,7) %>% rename(panos.id=sampleA,LogFC=logCPM,FDR=Geneid,PValue=FDR) %>% 
  mutate(ha_geneid=panos.rename(panos.id),Larval_Stage='L3') %>% filter(FDR<1e-4) %>% filter(ha_geneid!='LOC')
l4.gut=fread('./inputs/DE_tables/L4MG__C4__L4-UP.csv') %>% select(1,4,6,7) %>% rename(panos.id=sampleA,LogFC=logCPM,FDR=Geneid,PValue=FDR) %>% 
  mutate(ha_geneid=panos.rename(panos.id),Larval_Stage='L4') %>% filter(FDR<1e-4) %>% filter(ha_geneid!='LOC')                                                                                           

####################################### Make venn diagram ####################################### 
larval_stage_carcass_upregulated=rbindlist(list(l2.gut,l3.gut,l4.gut)) %>% merge(select(prot.key,-Length),by='ha_geneid') %>% select(ha_geneid,ha_prot_acc,Larval_Stage,LogFC,PValue,FDR)
write_xlsx(larval_stage_carcass_upregulated,'./outputs/Gut_DEGs/Larval_Stage_Carcass_Upregulated.xlsx')

v.l=list(l2.gut$ha_geneid,l3.gut$ha_geneid,l4.gut$ha_geneid)
names(v.l)=c('L2 Midgut \nEnriched','L3 Midgut \nEnriched','L4 Midgut \nEnriched')


life.ven=venn.diagram(v.l,filename=NULL,main='Life stage midgut specific genes \n \n',main.cex=6,cat.cex=2.5,cat.dist=.08,cex=3,fill=c('coral','cadetblue','gold'))#,cat.pos=c(315,60,180))
ggsave(life.ven,file='./outputs/Gut_DEGs/Life_stage_gut_specific_venn.tiff',device='tiff',width=14,height=12)
ggsave(life.ven,file='./outputs/Gut_DEGs/TableSX_Life_stage_gut_specific_venn.pdf',device='pdf',width=14,height=12)
for(i in list.files('./',full.names = T)){if(grepl('log',i)){file.remove(i)}}




common.genes=intersect(l2.gut$ha_geneid,l3.gut$ha_geneid) %>% intersect(l4.gut$ha_geneid) %>% writeLines('./outputs/gene_groups/Life_stage_common_upregulated.txt')
l2.specific=setdiff(l2.gut$ha_geneid,c(l3.gut$ha_geneid,l4.gut$ha_geneid)) %>% writeLines('./outputs/gene_groups/Life_stage_MidgutUp_L2_specific.txt')
l3.specific=setdiff(l3.gut$ha_geneid,c(l2.gut$ha_geneid,l4.gut$ha_geneid)) %>% writeLines('./outputs/gene_groups/Life_stage_MidgutUp_L3_specific.txt')
l4.specific=setdiff(l4.gut$ha_geneid,c(l3.gut$ha_geneid,l4.gut$ha_geneid)) %>% writeLines('./outputs/gene_groups/Life_stage_MidgutUp_L4_specific.txt')




####################################### FuzzyClusters for Life stage ####################################### 
set.seed(12346)
fuzzy.matrix=trans.sum %>% select(ha_geneid,matches('L.+MG_Artificial_trans_mean'),-matches("L5")) %>%  
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
ggsave(plot=gp,filename='./outputs/fuzzyClustering/Life_Stage_Custer_size_step.tiff',device='tiff',height=8,width=14)



####################################################  Make Fuzzy Plot #################################################### 
c.8 <- mfuzz(gut.exp.s,c=4,m=m1)
#c.8 <- mfuzz(gut.exp.s,c=6,m=m1)
#c.8 <- mfuzz(gut.exp.s,c=8,m=m1)
par(mar=c(1,1,1,1))
tiff('./outputs/fuzzyClustering/LifeStage_fuzzy.tiff',width=600,height=600)
fuzzy.plot=mfuzz.plot2(gut.exp.s,cl=c.8,mfrow=c(2,3),
                       time.labels = c('L2','L3','L4'),
                       xlab='Life Stage',x11=F,centre=T,min.mem=.4)
print(fuzzy.plot)
dev.off()



####################################################  Extract Fuzzy Genes #################################################### 
## annotate best cluster data and get best member
cluster=data.table(best.clust=c.8[[3]],ha_geneid=names(c.8[[3]]))
membership=c.8[[4]] %>% data.table()

for(i in 1:length(colnames(membership))){colnames(membership)[i]=paste0('A',colnames(membership)[i])} ### make column numeric
membership=membership %>% rowwise() %>% mutate(maximum=max(A1,A2,A3,A4)) %>% data.table() 
membership=cbind(cluster,membership)
fin.membership=membership %>% filter(maximum>.6) %>% arrange(desc(best.clust)) %>% data.table() 
mem.tab=fin.membership$best.clust %>% table()


####################################################  Annotate Fuzzy Clusters #################################################### 
key=data.table(L4_specific=2055,L2_specific=383,stage_going_down=238,stage_L3_down=960) %>% t() ### for .6 and 4 clusters

ml=list()
for(i in 1:4){
  n=dim(fin.membership[best.clust==i])[1]
  ml[[i]]=data.table(num=i,name=names(key[,][key[,]>n-100 & key[,]<n+100]))
}
final.fuzzy.key=rbindlist(ml)

cl_name=c()
for(i in fin.membership$best.clust){cl_name=c(cl_name,final.fuzzy.key[num==i]$name)}
fin.membership$cluster_name=cl_name



### take a peak which pH genes come up
#ph=fread('./inputs/acid_base/HelArm_pH_genes.csv') %>% filter(Category!='SLC_2') %>% mutate(ha_geneid=paste0('LOC',ha_geneid)) %>%
#merge(fin.membership,by='ha_geneid') %>% arrange(cluster_name)

####################################################  Write Fuzzy Gene groups #################################################### 

for(i in unique(fin.membership$cluster_name)){
  fil=fin.membership %>% filter(cluster_name==i)
  #writeLines(as.character(fil$ha_geneid),paste0('./outputs/gene_groups/Cluster_',as.character(i),'_Gene_List.txt'))
}


