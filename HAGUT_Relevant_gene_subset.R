######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(scales))
shhh(library(xlsx))


######################################## import data ######################################## 
setwd('/home/shanedenecke/Dropbox/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/Gene_Subset_expression')
merged.omics=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')

key=fread('./inputs/keys/Ha_master_key.csv') %>% select(Ha_OGS,ha_geneid) %>% mutate(ha_geneid=paste0('LOC',ha_geneid))
p450s=fread('./inputs/keys/Pearce_2017_Ha_Hz_equiv.csv') %>% select(`HaOGS2-id`,`function`) %>% rename(Ha_OGS=`HaOGS2-id`,gene_name=`function`) 



########## Table S6

trans.exp=trans.raw %>% select(-uni_id) %>% merge(merged.omics,by='ha_geneid')
write.xlsx(trans.exp,'./outputs/Gene_Subset_expression/TableS6_Transporter_expression_subset.xlsx',row.names = F)

merged.omics[ha_geneid %in% a$gene_name]

########## Table S7. Proteome upside
file.remove('./outputs/Gene_Subset_expression/TableS7_Proteome_spatial_compare.xlsx')
for(i in names(final.list)){
  sub=merged.omics[ha_geneid %in% final.list[[i]]]
  if(file.exists('./outputs/Gene_Subset_expression/TableS7_Proteome_spatial_compare.xlsx')){
    write.xlsx(sub,'./outputs/Gene_Subset_expression/TableS7_Proteome_spatial_compare.xlsx',append=T,sheetName = i,row.names = F)
  }else{
    write.xlsx(sub,'./outputs/Gene_Subset_expression/TableS7_Proteome_spatial_compare.xlsx',sheetName = i,row.names = F)
  }
  
}

########## Table S8

ph.merge=merge(ph,merged.omics,by='ha_geneid') %>% filter(Category %in% c('SLC_9','vATPase')) %>% select(-ha_prot_acc,-annot)

write.xlsx(ph.merge,'./outputs/Gene_Subset_expression/TableS9_pH_Gene_expression.xlsx')



################ Table S4 Proteome confusion matrix

a=merged.omics[ha_geneid %in% high.prot] %>% mutate(Category='High Proteome Low Transcriptome') 
b=merged.omics[ha_geneid %in% high.trans] %>% mutate(Category='High Transcriptome Low Proteome')

c=rbind(a,b)
fwrite(c,'./outputs/Gene_Subset_expression/TableS4_Transcriptome_Proteome_confusion.csv')
