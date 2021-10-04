#!/usr/bin/env Rscript

######################################## Import Packages ########################################
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
library(factoextra)

######################################## Set working directory ########################################
setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/proteomics_clean/',showWarnings = F)
######################################## Import proteome files ########################################
gel.free=fread('./inputs/proteomics/Ha_proteomics_L2-L5.csv') %>%
   select(Ha_prot,L3_prot,PM_prot,L4_prot,L2_prot) %>% filter_at(vars(contains('L')),any_vars(.=='Present'))

# key=fread('./inputs/keys/protein2gids.map',col.names = c('ha_prot_acc','ha_geneid','internal_id')) %>%
#   select(ha_geneid,ha_prot_acc) 
ncbi.key=fread('./inputs/keys/HelArm_NCBI_table.csv') %>% select(GeneID,`Protein product`,Length,`Protein Name`) %>% 
   rename(ha_geneid=GeneID,ha_prot_acc=`Protein product`,annot=`Protein Name`) %>% mutate(ha_geneid=paste0("LOC",ha_geneid)) %>%
   select(ha_geneid,ha_prot_acc) %>% unique.data.frame()


###################################### Functions ######################################################

rep.mean=function(x){
  sample=select(gel.based.clean,matches(x))
  return(rowMeans(sample,na.rm = T))
}

stderr=function(x){sd(x,na.rm = T)/sqrt(length(x))}

rep.std=function(x){
  sample=select(gel.based.clean,matches(x))
  return(apply(sample,1,stderr))
}


rep.sum=function(x){
  sample=select(gel.free.clean,matches(x))
  sample=apply(sample,2,as.numeric)
  return(rowSums(sample))
}




########################################  Process Raw proteome files ########################################  

########################################  Gel free ########################################  
gf.list=list()
#i="./inputs/proteomics/gel_free_csv/L1_3_mem.csv"
for(i in list.files('./inputs/proteomics/gel_free_csv/',full.names = T)){
  b=basename(i)
  if(grepl('PM_',b)){next}
  stage=gsub('(^L[0-9])_.+$','\\1',b)
  frac=gsub('^L[0-9].+_([a-z]{3}).csv','\\1',b)
  frac=ifelse(frac=='mem',yes='Membrane',no='Soluble')
  rep=gsub('^L[0-9].+([0-9])_.+.csv','\\1',b)
  gf.list[[i]]=fread(i) %>% filter(grepl('XP',Accession)) %>% select(Accession,`Score Mascot`) %>% mutate(Stage=stage,Fraction=frac,Replicate=rep)
}
gel.free.raw=rbindlist(gf.list) %>% rename(ha_prot_acc=Accession) %>% group_by(ha_prot_acc,Stage,Replicate) %>% summarize(Mascot_Score=sum(as.numeric(`Score Mascot`)))

gel.free.clean=gel.free.raw %>% #rename(Mascot_Score=`Score Mascot`) %>% 
  mutate(Sample=paste(Stage,Replicate,sep='_')) %>% select(-Stage,-Replicate) %>%
  pivot_wider(id_cols=ha_prot_acc,names_from=Sample,values_from=Mascot_Score) %>% data.table()
gel.free.clean[is.na(gel.free.clean)]=0

#fwrite(gel.free.clean,'./outputs/proteomics_clean/Gel_Free_Proteomics_Replicates.csv')



#### Calculate sums of peptides
samples=c('L2','L3','L4')
sam.sums=lapply(samples,rep.sum)

names(sam.sums)=samples

gel.free.sum=Reduce(cbind,sam.sums) %>% data.table()
colnames(gel.free.sum)=sapply(samples,paste0,'_prot_sum')

gel.free.sum$ha_prot_acc=gel.free.clean$ha_prot_acc

#fwrite(gel.free.sum,'./outputs/proteomics_clean/Gel_Free_Proteomics_Summary.csv')



###########################################   Gel based ########################################  
gb.list=list()
i="./inputs/proteomics/L5/L5_AM.f95.csv"
for(i in list.files('./inputs/proteomics/L5/',full.names = T)){
  b=basename(i)
  stage=gsub('(^L[0-9])_.+$','\\1',b)
  section=gsub('L[0-9]_([A-z]+).+$','\\1',b) %>% paste0(.,'G')
  #frac=gsub('^L[0-9].+_([a-z]{3}).csv','\\1',b)
  gb.list[[i]]=fread(i) %>% 
    rename_with(~gsub('nsaf_','',.x) %>% 
                  gsub('s_','Soluble_',.) %>% 
                  gsub('p_','Membrane_',.)) %>%
                  mutate(Stage=stage,
                         Section=section)
}
gel.based.raw=rbindlist(gb.list) %>% rename(ha_prot_acc=Protein_ID)

gel.based.clean=pivot_longer(gel.based.raw,cols=Soluble_1:Membrane_3,names_to="Fraction",values_to='Abundance') %>%
  mutate(Sample=paste0(Section,'_',Fraction)) %>% select(-Section,-Fraction) %>% 
  separate(Sample,into=c('Tissue','Fraction','Replicate'),sep='_') %>% 
  group_by(ha_prot_acc,Stage,Tissue,Replicate) %>% summarize(Abundance=sum(Abundance)) %>%
  mutate(Sample=paste(Tissue,Replicate,sep='_')) %>% select(-Stage,-Tissue,-Replicate) %>% 
  pivot_wider(id_cols=c(ha_prot_acc),names_from=Sample,values_from=Abundance) %>% 
  select(ha_prot_acc ,contains('FG'),contains('AMG'),contains('MMG'),contains('PMG'),contains('HG'))  %>% data.table()
gel.based.clean[is.na(gel.based.clean)]=0


#fwrite(gel.based.clean,'./outputs/proteomics_clean/Gel_Based_Proteomics_Replicates.csv')


samples=c('FG','AMG','MMG','PMG','HG')
sam.means=lapply(samples,rep.mean)
sam.std=lapply(samples,rep.std)

names(sam.means)=samples
names(sam.std)=samples

prot.summary1=Reduce(cbind,sam.means) %>% data.table()
colnames(prot.summary1)=sapply(samples,function(x) paste0('L5_',x,'_prot_mean'))

prot.summary2=Reduce(cbind,sam.std) %>% data.table()
colnames(prot.summary2)=sapply(samples,function(x) paste0('L5_',x,'_prot_std'))


gel.based.sum=cbind(prot.summary1,prot.summary2)
gel.based.sum$ha_prot_acc=gel.based.clean$ha_prot_acc

#fwrite(gel.based.sum,'./outputs/proteomics_clean/Gel_Based_Proteomics_Summary.csv')



#### Merge all proteomics


###### Summary

prot.sum=merge(gel.based.sum,gel.free.sum,by='ha_prot_acc',all=T) %>% merge(ncbi.key,by='ha_prot_acc',all=T)
prot.sum[is.na(prot.sum)]=0

dups2=prot.sum[duplicated(ha_geneid)]$ha_geneid

dups.table=prot.sum[ha_geneid %in% dups2]
nondups.table=prot.sum[!(ha_geneid %in% dups2)] %>% select(-ha_prot_acc)
cond=dups.table %>% select(-ha_prot_acc) %>% group_by(ha_geneid) %>% summarize_all(.funs=sum) %>% data.table()

prot.sum2=rbind(cond,nondups.table)

fwrite(prot.sum2,'./outputs/proteomics_clean/Final_Proteomics_Quant_Summary.csv')

### Presence Absence
prot.pres=prot.sum2 %>% select(ha_geneid,contains('mean'),contains('sum')) %>% rename_with(.fn=function(x) gsub('_mean','',x) %>% gsub('_sum','',.)) 
prot.pres[,2:ncol(prot.pres)][prot.pres[,2:ncol(prot.pres)]>0]="Present"
prot.pres[prot.pres=="0"]="Absent"

fwrite(prot.pres,'./outputs/proteomics_clean/Final_Proteomics_Presence.csv')



#### Repliacts
prot.rep=merge(gel.based.clean,gel.free.clean,by='ha_prot_acc',all=T) %>% merge(ncbi.key,by='ha_prot_acc',all=T) 
prot.rep[is.na(prot.rep)]=0
prot.rep=select(prot.rep,ha_geneid,ha_prot_acc,everything())
prot.rep[,3:ncol(prot.rep)]=lapply(prot.rep[,3:ncol(prot.rep)],as.numeric) %>% as.data.table()

dups2=prot.rep[duplicated(ha_geneid)]$ha_geneid

dups.table=prot.rep[ha_geneid %in% dups2]
nondups.table=prot.rep[!(ha_geneid %in% dups2)] %>% select(-ha_prot_acc)
cond=dups.table %>% select(-ha_prot_acc) %>% group_by(ha_geneid) %>% summarize_all(.funs=sum) %>% data.table()

prot.rep2=rbind(cond,nondups.table)


fwrite(prot.rep2,'./outputs/proteomics_clean/Final_Proteomics_Quant_Replicates.csv')
















pca.input=gel.based.clean %>% select(-ha_prot_acc) %>% as.matrix() %>% t()

#pca.input= select(gel.based.raw,Soluble_1:Membrane_3) %>% as.matrix() %>% t()
grouping=c(rep('FG',3),rep('AMG',3),rep('MMG',3),rep('PMG',3),rep('HG',3))
#grouping=c(rep('FG_sol',3),rep('FG_mem',3),rep('AMG_sol',3),rep('AMG_mem',3),rep('MMG_sol',3),rep('MMG_mem',3),rep('PMG_sol',3),rep('PMG_mem',3),rep('HG_sol',3),rep('HG_mem',3))
pal=c("red",'green','black','purple','yellow4',"blue",'orange','cyan','grey50','black','gold')
#pal=c('red1','red4','blue1','blue4','yellow1','yellow4','grey50','grey20','orange1','orange4')


PCA1=prcomp(pca.input)
pca=fviz_pca_ind(PCA1,
                 col.ind = grouping, 
                 geom='point',
                 palette = pal,
                 repel = TRUE,
                 addEllipses = TRUE, 
                 ellipse.type = "confidence",
                 legend.title = "Sample Groups",
                 title='L5 Section Proteomics PCA')
pca=pca+theme_bw()
pca=pca+theme(
  panel.grid=element_blank(),
  rect=element_rect(fill=NULL),
  text=element_text(size=18),
  axis.title=element_text(size=20),
  legend.title=element_text(size=24,face='bold'),
  plot.title = element_text(size=32,face='bold',hjust = 0.5))
print(pca)

ggsave(pca,filename='./outputs/PCA_plots/Spatial_proteomics_PCA.png',device='png',width=14,height=10)






pca.input.stage=prot.rep2 %>% select(contains('L2'),contains('L3'),contains('L4')) %>% 
  filter_at(vars(contains('L')),any_vars(.>0)) %>% as.matrix() %>% t()

grouping=c(rep('L2',3),rep('L3',3),rep('L4',3))
pal=c("red",'green','black','purple','yellow4',"blue",'orange','cyan','grey50','black','gold')


PCA1=prcomp(pca.input.stage)
pca=fviz_pca_ind(PCA1,
                 col.ind = grouping, 
                 geom='point',
                 palette = pal,
                 repel = TRUE,
                 addEllipses = TRUE, 
                 ellipse.type = "confidence",
                 legend.title = "Sample Groups",
                 title='Life Stage PCA')
pca=pca+theme_bw()
pca=pca+theme(
  panel.grid=element_blank(),
  rect=element_rect(fill=NULL),
  text=element_text(size=18),
  axis.title=element_text(size=20),
  legend.title=element_text(size=24,face='bold'),
  plot.title = element_text(size=32,face='bold',hjust = 0.5))
print(pca)

ggsave(pca,filename='./outputs/PCA_plots/Life_Stage_proteomics_PCA.png',device='png',width=14,height=10)












# ######################################## Import ########################################
# n0=list.files('./inputs/proteomics/L5/',full.names=T)
# n1=gsub(".f95.csv","",n0)
# n2=paste0(n1,'G')
# names(n2)=n0
# 
# 
# ########################################  Combine L5 data ########################################
# l=list()
# for(i in n0){
#   ind=fread(i)$Protein_ID
#   nam=n2[names(n2)==i]
#   df=data.table(Ha_prot=ind,junk='Present') %>% unique.data.frame()
#   colnames(df)=c('Ha_prot',nam)
#   l[[n2[names(n2)==i]]]=df
# }
# 
# final=Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Ha_prot", all = TRUE),l) ### merge all dataframes in list
# final[is.na(final)]='Absent'
# colnames(final)=basename(colnames(final))
# 
# ########################################  Merge with L2-L4 Data ########################################
# 
# total=merge(gel.free,final,by='Ha_prot',all=T) %>% rename(ha_prot_acc=Ha_prot) %>% merge(ncbi.key,by='ha_prot_acc',all=T) %>%
#   select(ha_geneid,ha_prot_acc,L2_prot,L3_prot,L4_prot,L5_FG,L5_AMG,L5_MMG,L5_PMG,L5_HG) %>% select(-ha_prot_acc) %>%
#   mutate(ha_geneid=paste0('LOC',ha_geneid))
# total$ha_geneid=gsub('LOC99','99',total$ha_geneid)
# total[is.na(total)]='Absent'
# 
# 
# ########################################  MERGE DUPLICATES ########################################
# 
# dups=total[ha_geneid %in% total[duplicated(ha_geneid)]$ha_geneid]
# total=total[!(ha_geneid %in% dups$ha_geneid)]
# l=list()
# for(i in unique(dups$ha_geneid)){
#   sub=dups[ha_geneid==i]
#   subprot=sub %>% select(-ha_geneid)
#   calls=apply(subprot,2, function(x) ifelse(('Present' %in% x),'Present','Absent'))
#   l[[i]]=cbind(select(sub,ha_geneid)[1],data.frame(as.list(calls)))
# }
# nodups=rbindlist(l)
# 
# final.final=rbind(total,nodups)
# colnames(final.final)=gsub('G$','G_prot',colnames(final.final))
# 
# fwrite(final.final,'./outputs/proteomics_clean/Ha_proteomics_L2-L5_Presence_Absence.csv')
# 
# 
