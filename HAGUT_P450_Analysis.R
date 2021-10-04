#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(scales))
shhh(library(ggsci))


######################################## import data ######################################## 
#setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
setwd('C:/Users/shane/Dropbox/Shane/omics_projects/Ha_midgut_atlas/')
dir.create('./outputs/Gene_Subset_expression/',showWarnings = F)
######################################## Clean data ######################################## 
merged.omics=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv',header=T,skip = 1) %>% rename('ha_geneid'=`NCBI Gene ID`,ha_prot_acc=`NCBI Protein ID (Longest Isoform)`)
key=fread('./inputs/keys/Ha_master_key.csv') %>% select(Ha_OGS,ha_geneid) 
#merged.omics=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')



################################# ALLLLL DETOX GENES ##############

special=fread('./inputs/keys/Pearce_2017_Ha_Hz_equiv.csv') %>% select(`HaOGS2-id`,`function`) %>% rename(Ha_OGS=`HaOGS2-id`,gene_name=`function`) %>% 
  filter(grepl('CYP|ABC|CCE|GST|Tryp|Chym',gene_name)) %>% merge(key,by='Ha_OGS',all.x=T) %>% unique.data.frame() %>% 
  filter(complete.cases(.)) %>% mutate(ha_geneid=paste0('LOC',ha_geneid)) %>%
  merge(merged.omics,by='ha_geneid')  %>% select(gene_name,matches('mean|std')) %>% select(-contains('prot'))
 


avg=melt(select(special,gene_name,matches('mean')),id.vars='gene_name',variable.name='Tissue',value.name='Expression') #%>% mutate(Expression=Expression+1)
std=melt(select(special,gene_name,matches('std')),id.vars='gene_name',variable.name='Tissue',value.name='std_err')

p2=cbind(avg,select(std,std_err)) %>%
  mutate(Diet=ifelse(grepl('Plant',Tissue),'Plant','Artificial')) %>% 
  mutate(Stage=gsub('(^L[0-9])_.+$','\\1',Tissue)) %>% 
  mutate(Tissue=gsub('L[0-9]_(.+G)_[Artificial|Plant].+$','\\1',Tissue)) %>%
  mutate(family=gsub('Harm(.{3}).+$','\\1',gene_name)) %>%
  filter(nchar(family)<4) %>% 
  mutate(family=gsub('Try','Trypsin',family)) %>%
  mutate(family=gsub('Chy','Chymotrypsin',family)) %>%
  filter(nchar(Tissue)<4)


p2$Tissue=factor(p2$Tissue,levels=c('FG','AMG','MMG','PMG','HG'))
p2$gene_name=gsub('Harm','',p2$gene_name)
p2=p2[complete.cases(p2)]

######################################## Subset out Cyp6 plant fed samples ######################################## 
plot.data=p2 %>% 
  filter((grepl('ABCB|ABCC|ABCG|CYP6|CYP9|GSTD|GSTT|CCE|Tryp|Chym',gene_name))) %>%
  filter(Stage=='L5') %>% 
  filter(Tissue!='C')
#  filter(grepl('HarmCYP6|HarmCYP9',gene_name)) %>% 
  



#################################### RAW #################################################
gp=ggplot(plot.data,aes(x=Tissue,y=Expression,fill=Diet))
gp=gp+ stat_boxplot(geom ='errorbar', width = 0.6)
gp=gp+geom_boxplot(width=.6)
#gp=gp+geom_dotplot(binaxis='y',stackdir='center')#,aes(fill=gene_name))
gp=gp+facet_wrap(facets=vars(family),nrow=2)
gp=gp+scale_fill_manual(values=c('grey50','yellow4','blue','grey50','orange',1:30))
gp=gp+labs(y='Expression TPM')
gp=gp+ scale_y_continuous(trans = pseudo_log_trans(),breaks = c(2,5,10,20,40,100,200,400))
gp=gp+ggtitle('Plant vrs. Diet Comparisons')
gp=gp+theme_bw()
gp=gp+theme(legend.title=element_text(size=24,face='bold'),
            legend.text=element_text(size=20,face='bold'),
            axis.text=element_text(size=20),axis.title=element_text(size=24),
            #legend.position='none',
            strip.text = element_text(size=28,hjust = 0.5),
            plot.title = element_text(size=36,hjust = 0.5))

print(gp)
ggsave(filename='./outputs/Gene_Subset_expression/Special_Diet_plant_raw',gp,device='png',width=12,height=8)



#################################### RATIO PLOT #################################################
ratio.data=spread(select(plot.data,-std_err),key=Diet,value=Expression) %>% mutate(plant_ratio=log2(Plant/Artificial))
gp=ggplot(ratio.data,aes(x=Tissue,y=plant_ratio,color=Tissue,shape=Tissue))
#gp=gp+ stat_boxplot(geom ='errorbar', width = 0.6,color='black')
#gp=gp+geom_boxplot(width=.6,fill='grey',color='black',outlier.shape=NA)#gp=gp+geom_dotplot(binaxis='y',stackdir='center')
gp=gp+stat_summary(geom='crossbar',fun='mean',size=.5,width=.4)
gp=gp+stat_summary(geom='errorbar',fun.data='mean_cl_normal',size=1.5,width=.6)
gp=gp+geom_point(position='jitter',alpha=.5,size=2)
#gp=gp+geom_violin()
gp=gp+scale_color_aaas()
gp=gp+geom_hline(yintercept = 0,linetype=2)
gp=gp+facet_wrap(facets=vars(family),nrow=2)
gp=gp+labs(y='Artificial <--  log2(Fold Change)   --> Plant\n')
gp=gp+ggtitle('Plant vrs. Diet Comparisons')
gp=gp+theme_bw()
gp=gp+theme(legend.title=element_text(size=24,face='bold'),
            legend.text=element_text(size=20,face='bold'),
            axis.text=element_text(size=20),axis.title=element_text(size=24),
            panel.grid = element_blank(),
            legend.position='none',
            strip.text = element_text(size=28,hjust = 0.5),
            plot.title = element_text(size=36,hjust = 0.5))

print(gp)
ggsave(filename='./outputs/Gene_Subset_expression/Special_Diet_Plant_ratio.pdf',gp,device='pdf',width=14,height=8)

l=list()
for(i in unique(ratio.data$family)){
  for(j in unique(ratio.data$Tissue)){
    sub=ratio.data[family==i & Tissue==j]
    test=sub$plant_ratio[!is.infinite(sub$plant_ratio) & !is.na(sub$plant_ratio)]
    l[[paste(i,j,sep='__')]]=t.test(test)$p.val
  }}
  


########################## Tables
s.10=fread('./inputs/keys/Pearce_2017_Ha_Hz_equiv.csv') %>% rename(Ha_OGS=`HaOGS2-id`,gene_name=`function`) %>% select(Ha_OGS,gene_name) %>%
  filter(grepl('CYP|GST|CCE|ABC|Try|Chy',gene_name)) %>% 
  filter(nchar(gene_name)<20) %>% merge(key,by='Ha_OGS') %>%
  unique.data.frame() %>% 
  mutate(ha_geneid=paste0('LOC',ha_geneid)) %>%
  merge(merged.omics,by='ha_geneid') %>%
  mutate(family=gsub('Harm(.{3}).+$','\\1',gene_name)) %>%
  mutate(family=gsub('Try','TRYPSIN',family)) %>%
  mutate(family=gsub('Chy','CHYMOTRYPSIN',family))
  
fwrite(s.10,'./outputs/Gene_Subset_expression/TableS10_Detox_Trypsin.csv')
