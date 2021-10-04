#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(factoextra))
library(gridExtra)

######################################## Directories ######################################## 
setwd('C:/Users/shane/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
#setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/omics_clean/',showWarnings = F)
dir.create('./outputs/PCA_plots',showWarnings = F)

######################################## Import raw files ######################################## 
ha.gid.key=fread('./inputs/keys/protein2gids.map')
colnames(ha.gid.key)=c('ha_prot_acc','ha_geneid','Geneid')
exclude=c('LOC110378454','LOC110379977')
ncbi.key=fread('./inputs/keys/HelArm_NCBI_table.csv') %>% select(GeneID,`Protein product`,Length,`Protein Name`) %>% 
  rename(ha_geneid=GeneID,ha_prot_acc=`Protein product`,annot=`Protein Name`)
raw.counts.mg=fread('./inputs/transcriptome/Ha_Midgut_Counts.txt') 
#proteomics=fread('./outputs/proteomics_clean/Ha_proteomics_L2-L5_Jan_2020.csv') #%>% select(-ha_geneid)
proteomics.sum=fread('./outputs/proteomics_clean/Final_Proteomics_Quant_Summary.csv')
proteomics.rep=fread('./outputs/proteomics_clean/Final_Proteomics_Quant_Replicates.csv')
proteomics.pres=fread('./outputs/proteomics_clean/Final_Proteomics_Presence.csv')
########################################Import functions ######################################## 


col.clean=function(x){
  x1=gsub('d_','G_Artificial_',x)
  x2=gsub('p_','G_Plant_',x1)
  x3=gsub('(L[2-4])_([0-9]$)','\\1_MG_Artificial_\\2',x2)
  x4=gsub('(L[2-4]_C)','\\1_Artificial',x3)
  x5=gsub('L5_C_','L5_C_Artificial_',x4)
  x6=gsub('L1_','L1_Artificial_',x5)
  return(x6)
}


stderr=function(x){sd(x)/sqrt(length(x))}

rep.mean=function(x){
  sample=select(tpm.rep,matches(x))
  return(rowMeans(sample,na.rm = T))
}

stderr=function(x){sd(x,na.rm = T)/sqrt(length(x))}

rep.std=function(x){
  sample=select(tpm.rep,matches(x))
  return(apply(sample,1,stderr))
}

######################################## Get Unigene set of 1 protein sequence per proteome 

uniprot=ncbi.key %>% group_by(ha_geneid) %>% summarize(ha_prot_acc=ha_prot_acc[Length==max(Length)][1],
                                                       Length=Length[Length==max(Length)][1],
                                                       annot=annot[Length==max(Length)][1]) %>% data.table() %>% 
  .[['ha_prot_acc']]

######################################## Process TPM Data by sample ######################################## 

### clean counts file
raw.counts=merge(raw.counts.mg,ha.gid.key,by='Geneid') %>% 
  merge(ncbi.key,by=c('ha_prot_acc','ha_geneid'),all=T) %>%
        filter(ha_prot_acc %in% uniprot) %>% 
  select(ha_prot_acc:ha_geneid,Length,annot,everything(),-Geneid,-contains('L5_C_'),-contains('L1')) %>% data.table()
raw.counts$ha_geneid=paste0("LOC",raw.counts$ha_geneid)
raw.counts=raw.counts %>% filter(!(ha_geneid %in% exclude))



### Counts to TPM
cols=colnames(raw.counts)  
iter=grep("L[0-9]",cols) 
colnames(raw.counts)=col.clean(colnames(raw.counts))
colnames(raw.counts)=gsub('(L.+$)','\\1_Counts',colnames(raw.counts))
normlen=raw.counts$Length/1000 

full.transcriptomes=raw.counts
for(i in iter){
  relevant.column=raw.counts[[i]]
  rpk=relevant.column/normlen
  sc=sum(rpk,na.rm = T)/1000000
  tpm=rpk/sc
  new.col=paste0(col.clean(cols[i]),'_TPM')
  full.transcriptomes[[new.col]]=tpm
}

#### Get full transcriptomes all reps
full.trans=full.transcriptomes %>% merge(proteomics.pres,by='ha_geneid')
fwrite(full.trans,'./outputs/omics_clean/Ha_Transcriptomics_Proteomes_full_replicates.csv')


tpm.rep=full.transcriptomes %>% select(ha_geneid,ha_prot_acc,contains('TPM')) 
colnames(tpm.rep)=gsub('_TPM','',colnames(tpm.rep))
#fwrite(tpm.rep,'./outputs/omics_clean/Ha_Transcriptomics_Replicate.csv')

#full.raw=select(raw.counts,-ha_prot_acc,-Length,-annot) %>% merge(tpm.data,by='ha_geneid')
#fwrite(full.raw,'./outputs/omics_clean/TableSY_FullRAW_data_by_sample.csv')

######################################## Melt data for ggplot2 graphs ######################################## 

#m=melt(tpm.rep,id.vars=c('ha_prot_acc','ha_geneid','Length','annot'),variable.name ='Tissue',value.name='TPM') %>% 
m=melt(tpm.rep,id.vars=c('ha_geneid'),variable.name ='Tissue',value.name='TPM') %>% 
  mutate(Diet_type=ifelse(grepl('Plant',Tissue),"Plant",'Artificial')) %>% data.table()
m$Replicate=gsub('^.+_([0-9])','\\1',m$Tissue)
m$Tissue=gsub('(^.+)_Artificial_[0-9]','\\1',m$Tissue)
m$Tissue=gsub('(^.+)_Plant_[0-9]','\\1',m$Tissue)
fwrite(m,'./outputs/omics_clean/Ha_Transcripomics_TPM_Replicate_Melt.csv')



######################################## Summarize Data by Sample ######################################## 
samples=c(
          'L2_MG_Artificial_','L2_C_Artificial_','L3_MG_Artificial_','L3_C_Artificial_','L4_MG_Artificial_','L4_C_Artificial_',
          'L5_FG_Artificial_','L5_FG_Plant_','L5_AMG_Artificial_','L5_AMG_Plant_',
          'L5_MMG_Artificial_','L5_MMG_Plant_', 'L5_PMG_Artificial_','L5_PMG_Plant_', 'L5_HG_Artificial_','L5_HG_Plant_')

sam.means=lapply(samples,rep.mean)
sam.std=lapply(samples,rep.std)

names(sam.means)=samples
names(sam.std)=samples

tpm.summary1=Reduce(cbind,sam.means) %>% data.table()
colnames(tpm.summary1)=sapply(samples,paste0,'trans_mean')

tpm.summary2=Reduce(cbind,sam.std) %>% data.table()
colnames(tpm.summary2)=sapply(samples,paste0,'trans_std')

tpm.summary=cbind(tpm.summary1,tpm.summary2)

tpm.summary$ha_geneid=raw.counts$ha_geneid
#tpm.summary$ha_prot_acc=raw.counts$ha_prot_acc
tpm.summary$annot=raw.counts$annot

#colnames(tpm.summary)=gsub('([A-Z])mean','\\1_mean',colnames(tpm.summary))

tpm.summary=tpm.summary  %>% select(ha_geneid,annot,everything())


######################################## Write data ######################################## 



merged.omics.sum=merge(tpm.summary,proteomics.pres,by='ha_geneid',all=T) %>% 
  merge(select(full.transcriptomes,ha_prot_acc,ha_geneid),by='ha_geneid') %>%
  select(ha_geneid,ha_prot_acc,everything()) #%>% rename(L5_C_Artificial_std=L5_C_Artificialstd)
merged.omics.sum[is.na(merged.omics.sum)]=0
fwrite(merged.omics.sum,'./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')

# merged.omics.rep=merge(tpm.rep,proteomics.rep,by='ha_geneid',all=T) %>% select(ha_geneid,everything()) #%>% rename(L5_C_Artificial_std=L5_C_Artificialstd)
# merged.omics.rep[is.na(merged.omics.rep)]=0
# fwrite(merged.omics.rep,'./outputs/omics_clean/Ha_proteomics_transcriptomics_Replicates.csv')

# merged.omics.pres=merge(tpm.summary,proteomics.pres,by='ha_geneid',all=T) %>% select(ha_geneid,everything()) #%>% rename(L5_C_Artificial_std=L5_C_Artificialstd)
# merged.omics.pres[is.na(merged.omics.pres)]=0
# fwrite(merged.omics.pres,'./outputs/omics_clean/Ha_proteomics_transcriptomics_Presence.csv')


######################## Create PCA Plots

### PCA plots function
pca.plot=function(data,grouping,title,pal=c("red","blue",'green','black','purple','yellow4','orange')){
  data=
  
  PCA1=prcomp(data)
  
  pca=fviz_pca_ind(PCA1,
                  col.ind = grouping, 
                  #geom='point',
                  palette = pal,
                  repel = TRUE,
                  addEllipses = TRUE, 
                  ellipse.type = "confidence",
                  legend.title = "Sample Groups",
                  title=gsub('_',' ',title))
  pca=pca+theme_bw()
  pca=pca+theme(
    panel.grid=element_blank(),
    rect=element_rect(fill=NULL),
    text=element_text(size=18),
    axis.title=element_text(size=20),
    legend.title=element_text(size=24,face='bold'),
    plot.title = element_text(size=32,face='bold',hjust = 0.5))
  return(pca)

  ggsave(filename =paste0('./outputs/PCA_plots/',title,'_PCA.tiff'),device='tiff',height=7,width=14)
}

tpm.pca=tpm.rep[complete.cases(tpm.rep)]



#### Life stage comparisons
stage.table=tpm.pca %>% select(L2_MG_Artificial_1:L4_C_Artificial_4) %>% rename_with(.fn=function(x) str_replace_all( x, "_Artificial_", "" )) %>% as.matrix() %>% t()
stage.groups=c(rep('L2 Midgut',4),rep('L2 Carcass',4),rep('L3 Midgut',4),rep('L3 Carcass',4),rep('L4 Midgut',4),rep('L4 Carcass',4))
pca.plot(data=stage.table,grouping=stage.groups,title='Larval_Stage_Comparison')

### Artificial diet comparisons
l5.art.table=tpm.pca %>% select(matches('L5_.+Artificial')) %>% rename_with(.fn=function(x) str_replace_all( x, "_Artificial_", "" )) %>% as.matrix() %>% t()
l5.art.groups=c(rep('Foregut',3),rep('Anterior Midgut',3),rep('Middle Midgut',3),rep('Posterior Midgut',3),rep('Hindgut',3))
pca.plot(data=l5.art.table,grouping=l5.art.groups,title='Artificial_Fed_L5_Gut_Compartments')


### Plant fed L5 samples samples
l5.plant.table=tpm.pca %>% select(matches('L5_.+Plant')) %>% rename_with(.fn=function(x) str_replace_all(x, "_Plant_", "" )) %>% as.matrix() %>% t()
l5.plant.groups=c(rep('Foregut',4),rep('Anterior Midgut',4),rep('Middle Midgut',4),rep('Posterior Midgut',4),rep('Hindgut',4))
pca.plot(data=l5.plant.table,grouping=l5.plant.groups,title='Plant_Fed_L5_Gut_Compartments',pal=c("red",'green','black','purple','yellow4',"blue",'orange'))


l5.plant.art.table=tpm.pca %>% select(matches('L5_')) %>% select(-contains("_C_")) %>% as.matrix() %>% t()
l5.plant.art.groups=c(rep('FG_Artificial',3),rep('AMG_Artificial',3),rep('MMG_Artificial',3),rep('PMG_Artificial',3),rep('HG_Artificial',3),
                  rep('FG_Plant',4),rep('AMG_Plant',4),rep('MMG_Plant',4),rep('PMG_Plant',4),rep('HG_Plant',4))
cols=c("red1",'red4','green1','green4','black','grey30','purple1','purple4','yellow1','yellow4')

pca.plot(data=l5.plant.art.table,grouping=l5.plant.art.groups,title='Plant_vs._Art_L5_Gut_Compartments',pal=cols)



l5.fg=tpm.pca %>% select(matches('FG')) %>% as.matrix() %>% t()
l5.fg.groups=c(rep('FG_Artificial',3),rep('FG_Plant',4))
cols=c("yellow4",'green4')
l5.fg.pca=pca.plot(data=l5.fg,grouping=l5.fg.groups,title='FG_Diet_Comparison',pal=cols)


l5.AMG=tpm.pca %>% select(matches('AMG')) %>% as.matrix() %>% t()
l5.AMG.groups=c(rep('AMG_Artificial',3),rep('AMG_Plant',4))
cols=c("yellow4",'green4')
l5.amg.pca=pca.plot(data=l5.AMG,grouping=l5.AMG.groups,title='AMG_Diet_Comparison',pal=cols)


l5.MMG=tpm.pca %>% select(matches('MMG')) %>% as.matrix() %>% t()
l5.MMG.groups=c(rep('MMG_Artificial',3),rep('MMG_Plant',4))
cols=c("yellow4",'green4')
l5.mmg.pca=pca.plot(data=l5.MMG,grouping=l5.MMG.groups,title='MMG_Diet_Comparison',pal=cols)

l5.PMG=tpm.pca %>% select(matches('PMG')) %>% as.matrix() %>% t()
l5.PMG.groups=c(rep('PMG_Artificial',3),rep('PMG_Plant',4))
cols=c("yellow4",'green4')
l5.pmg.pca=pca.plot(data=l5.PMG,grouping=l5.PMG.groups,title='PMG_Diet_Comparison',pal=cols)

l5.HG=tpm.pca %>% select(matches('HG')) %>% as.matrix() %>% t()
l5.HG.groups=c(rep('HG_Artificial',3),rep('HG_Plant',4))
cols=c("yellow4",'green4')
l5.hg.pca=pca.plot(data=l5.HG,grouping=l5.HG.groups,title='HG_Diet_Comparison',pal=cols)


plant.diet=grid.plot=grid.arrange(l5.fg.pca,
                       l5.amg.pca,
                       l5.mmg.pca,
                       l5.pmg.pca,
                       l5.hg.pca)

ggsave(plot=plant.diet,filename = './outputs/PCA_plots/Plant_Diet_Combined_PCA.pdf',width=15,height=12)


