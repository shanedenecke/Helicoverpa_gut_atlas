#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(scales))
shhh(library(RColorBrewer))
shhh(library(edgeR))



######################################## import data ######################################## 
#setwd('/home/shanedenecke/Dropbox/omics_projects/Ha_midgut_atlas')
setwd('C:/Users/shane/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/artificial_plant_compare',showWarnings = F)


######################################## Import raw files ######################################## 
ha.gid.key=fread('./inputs/keys/protein2gids.map')
colnames(ha.gid.key)=c('ha_prot_acc','ha_geneid','Geneid')
ncbi.key=fread('./inputs/keys/HelArm_NCBI_table.csv') %>% select(GeneID,`Protein product`,Length,`Protein Name`) %>% 
  rename(ha_geneid=GeneID,ha_prot_acc=`Protein product`,annot=`Protein Name`)
raw.counts.mg=fread('./inputs/transcriptome/Ha_Midgut_Counts.txt')
#proteomics=fread('./outputs/proteomics_clean/Ha_proteomics_L2-L5_Jan_2020.csv') #%>% select(-ha_geneid)



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

######################################## Get Unigene set of 1 protein sequence per proteome 

uniprot=ncbi.key %>% group_by(ha_geneid) %>% summarize(ha_prot_acc=ha_prot_acc[Length==max(Length)][1],
                                                       Length=Length[Length==max(Length)][1],
                                                       annot=annot[Length==max(Length)][1]) %>% data.table() %>% 
  .[['ha_prot_acc']]


######################################## clean counts file ######################################## 
raw.counts=merge(raw.counts.mg,ha.gid.key,by='Geneid') %>% 
  merge(ncbi.key,by=c('ha_prot_acc','ha_geneid')) %>%
  filter(ha_prot_acc %in% uniprot) %>% 
  select(ha_prot_acc:ha_geneid,Length,annot,everything(),-Geneid) %>% data.table()
raw.counts$ha_geneid=paste0("LOC",raw.counts$ha_geneid)

colnames(raw.counts)=col.clean(colnames(raw.counts))


########################################  Run Edge R ######################################## 
tissues=c('FG','AMG','MMG','PMG','HG')
counts.de=raw.counts %>% select(-Length,-annot,-ha_prot_acc)

de.list=list()
de.list2=list()
for(i in tissues){
  sub.de=select(counts.de,ha_geneid,contains(i))
  art.cols=colnames(sub.de)[grepl('Artificial',colnames(sub.de))]
  plant.cols=colnames(sub.de)[grepl('Plant',colnames(sub.de))]
  
  art.name=paste0(i,"_Artificial")
  plant.name=paste0(i,"_Plant")
  
  
  sub.de=select(sub.de,ha_geneid,all_of(art.cols),all_of(plant.cols))
    
  edge.de=sub.de
  rownames(edge.de)=raw.counts$ha_geneid
  edge.de$ha_geneid=NULL
  
  group=factor(c(rep(art.name,length(art.cols)),rep(plant.name,length(plant.cols))))
  delist=DGEList(counts=edge.de,group=group)
  
  ### filter
  keep=filterByExpr(delist)
  delist=delist[keep, , keep.lib.sizes=FALSE]
  
  
  ## Add normalization
  delist=calcNormFactors(delist)
  delist = estimateCommonDisp(delist)
  delist = estimateTagwiseDisp(delist)
  
  ##perform test
  et = exactTest(delist, pair=c(art.name, plant.name))
  tTags = topTags(et,n=NULL)
  
  output=tTags$table 
  row.index=rownames(output) %>% as.numeric()
  genes=raw.counts[row.index]$ha_geneid
  output$ha_geneid=genes
  
  signif=data.table(sampleA=art.name, sampleB=plant.name, output) %>% 
    merge(sub.de,by='ha_geneid',all=T) %>%
    arrange(FDR) %>%  
    filter(FDR<1e-3 & abs(logFC)>2) %>% mutate(comparison=i) %>%
    data.table() 
  
  de.list[[i]]=signif
  
  plantup=signif[logFC>2]$ha_geneid 
  artup=signif[logFC<2]$ha_geneid
  de.list2[[paste0(i,'art')]]=artup
  de.list2[[paste0(i,'plant')]]=plantup
  
  
  
  plantup%>% writeLines(paste0('./outputs/gene_groups/L5_',i,'_Plant_Enriched.txt'))
  artup %>% writeLines(paste0('./outputs/gene_groups/L5_',i,'_Artificial_Enriched.txt'))
  
}

art.plant.compare=rbindlist(de.list) %>% 
  merge(select(raw.counts,ha_geneid,annot),by='ha_geneid') %>% 
  arrange(comparison,FDR) %>%
  select(ha_geneid,annot,comparison,logFC,FDR)
fwrite(art.plant.compare,'./outputs/artificial_plant_compare/Artificial_plant_diet_compare_DE.csv')


########################### Split DE lists for GO term analysis  ########################### 
a=intersect(de.list2$HGplant,de.list2$FGplant) %>% intersect(de.list2$AMGplant) %>% intersect(de.list2$MMGplant) %>% intersect(de.list2$PMGplant)
