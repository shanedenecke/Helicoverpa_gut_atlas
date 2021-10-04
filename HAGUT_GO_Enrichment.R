#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh=suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))

######################################## Directories ######################################## 
#setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
setwd('C:/Users/shane/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
source('./scripts/HAGUT_Functions.R')
dir.create('outputs/GO',showWarnings = F)


######################################## Import raw files ######################################## 
ha.key.raw=data.table::fread('./inputs/keys/Ha_master_key.csv') %>% mutate(ha_geneid=paste0('LOC',ha_geneid)) 
ha.key=data.table::fread('./inputs/keys/Ha_master_key.csv') %>% mutate(ha_geneid=paste0('LOC',ha_geneid)) %>% select(ha_geneid,Ha_OGS) %>% unique.data.frame()

ha.omics.full=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')

ha.go=fread('./inputs/keys/Pearce_2017_GO.csv',header=F,col.names = c('Ha_OGS','code_name'))
go.key=fread('./inputs/keys/GO_key_clean.txt',col.names = c('code_name','GO_name','GO_cat'))
de.list=lapply(list.files('./outputs/gene_groups/',full.names = T),readLines)
names(de.list)=list.files('./outputs/gene_groups/') %>% gsub('_Gene_List.txt','',.)
go.full=merge(ha.go,go.key,by='code_name') %>% merge(ha.key,by='Ha_OGS') 


######################################## Functions ######################################## 

freq.table=function(x){
  #x=go.full
  total=x$ha_geneid %>% unique() %>% length()
  l=list()
  for(i in unique(x$code_name)){
    sub=x[code_name==i]
    num.genes=unique(sub$ha_geneid) %>% length()
    freq=num.genes/total
    l[[i]]=data.table(code_name=i,code_count=num.genes,frequency=freq,total=total)
  }
  return(rbindlist(l))
} 



vlepo.enrich=function(trial,omic){
  #trial=go.de.freq[[1]]
  #omic=goterms.freq
  if(nrow(trial)>0){
    
    fischer_p=c()
    test.list=list()
    ref.list=list()
    for(i in trial$code_name){
      test=trial[trial$code_name==i,c('code_count','total')] %>% as.numeric()
      ref=omic[omic$code_name==i,c('code_count','total')] %>% as.numeric()
      p.value=fisher.test(matrix(c(test[1],ref[1],test[2],ref[2]),nrow=2,ncol=2),alternative="greater")$p.value ## 3extract P value from fischer test
      fischer_p=c(fischer_p,p.value)
      
      names(test)=c('test_N','test_total')
      test.list[[i]]=as.data.frame(t(test))
      
      names(ref)=c('ref_N','ref_total')
      ref.list[[i]]=as.data.frame(t(ref))
    }
    
    test.add=rbindlist(test.list)
    ref.add=rbindlist(ref.list)
    trial$fischP=fischer_p 
    trial$fdr=p.adjust(trial$fischP,method='fdr')
    combined=cbind(test.add,ref.add,trial)
    enriched= combined %>% arrange(fdr) %>% filter(fischP<.05) %>% filter(test_N>5) %>% select(-frequency,-code_count) %>% 
      mutate(test_freq=test_N/test_total, ref_freq=ref_N/ref_total) %>% 
      select(code_name,fischP,fdr,everything()) %>% data.table()
    return(enriched)
  }
}


gene.subset=function(geneset,database){
  return(database[ha_geneid %in% geneset])
}


######################################## Perform GO enrichment ######################################## 


### get freuqnency table for entiere genome
goterms.freq=go.full %>% freq.table()

### extract GO terms and annotations for each gene 
de.annot=lapply(de.list,function(x) ha.omics.full[ha_geneid %in% x] %>% select(ha_geneid,annot)) ### get annotations for each gene in each de list
go.de=lapply(de.list, function(x) go.full[ha_geneid %in% x]) ### Extract GO terms from Pearce reference for each gene in each de list

### Get frequency table for each GO term in each DE list
go.de.freq=list() ### set up empty list
for(i in names(go.de)){go.de.freq[[i]]=freq.table(go.de[[i]])}

#### perform enrichment analysis
goterms.enrich=lapply(go.de.freq,vlepo.enrich,goterms.freq)

######################################## Annotate GO enrichment ######################################## 
l=list()
for(i in names(goterms.enrich)){
  enrich=goterms.enrich[[i]]
  if(is.null(enrich)){next} ### if there are no enriched go terms then skip
  annot=merge(enrich,go.key,by='code_name',all.y=F,sort=F)
  l[[i]]=data.table(annot[!duplicated(annot)],category=i)
}

######################################## Filter and Write GO enrichment ######################################## 
significant.go=rbindlist(l) %>% #### All GO terms 
  filter(fdr<5e-3)  %>% ### FDR threshold
  filter(GO_cat!='cellular_component')  %>% ### Set category if you want
  #select(code_name,fdr,GO_name,category) %>% 
  rename(`GO code`=code_name,`False Discovery Rate`=fdr,`GO Name`=GO_name,Comparison=category)


relevant.comparisons=significant.go %>% 
  filter(grepl('Cluster',Comparison) | grepl('MG_unique',Comparison) | grepl('Life_stage',Comparison))



fwrite(significant.go,'./outputs/GO/ALL_ENRICHED_GO.csv')


### This function is for manually checking which genes are present in a given group and GO term
 go_viz=function(group,code){
   ha.key.sub=select(ha.key.raw,ha_geneid,annot)
   temp=go.de[[group]][code_name==code]
   m2=merge(temp,ha.key.sub,by='ha_geneid') %>% .[!duplicated(ha_geneid)]
   return(m2)
 }
go_viz("Midgut_common_proteome.txt",'GO:0004252')



