#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(scales))
shhh(library(ggsci))

######################################## import data ######################################## 
setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
#setwd('C:/Users/shane/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/Trans_prot_compare',showWarnings = F)

merged.omics.pres=fread('./outputs/omics_clean/Ha_proteomics_transcriptomics_Presence.csv') %>% filter(complete.cases(.)) %>% 
  select(-matches('std'),-matches('Plant'))

merged.omics.sum=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv') %>% filter(complete.cases(.)) %>% 
  select(-matches('std'),-matches('Plant'))

     
prot.pos=merged.omics.pres %>% filter_at(vars(matches('prot$')),any_vars(.=='Present')) %>% .[['ha_geneid']]
merged.omics.pres=merged.omics.pres %>% mutate(overall_proteome=ifelse((ha_geneid %in% prot.pos),'Present','Absent'))



low.thresh=2
high.thresh=2000

######################################## Graphing function ########################################
omic.compare=function(x){
  #x='L2'
  sub=merged.omics.pres %>% select(ha_geneid,matches(x),-matches('_C_'))
  l=list()

  prot=colnames(sub)[grepl('prot',colnames(sub))]
  trans=colnames(sub)[grepl('mean',colnames(sub))]
  sub$trans_adjust=sub[[trans]]+1

  sub2=sub %>% `colnames<-`(c('ha_geneid','transcriptome','proteome','transcriptome_adjust')) %>% mutate(sample=x)
  l[[x]]=sub2

  ###extract outliers
  high.prot=sub %>% filter(get(prot)=='Present') %>% filter(get(trans)<low.thresh) %>% .[['ha_geneid']]
  high.trans=sub %>% filter(get(prot)=='Absent') %>% filter(get(trans)>high.thresh) %>% .[['ha_geneid']]

  #### Make plot
  gp=ggplot(data=sub,aes_string(x=prot,y='trans_adjust',fill='prot'))
  gp=gp+geom_violin()
  gp=gp+scale_fill_aaas()
  gp=gp+labs(x='\nProteome Call',y='Expression Level (TPM)\n')
  gp=gp+ggtitle(paste0(x,' Transcriptomic and Proteomic Comparison'))
  gp=gp+scale_y_continuous(trans = log2_trans(),breaks=c(2,5,10,20,50,100,200,400,800,1500,3000))
  gp=gp+theme_bw()
  gp=gp+theme(legend.position='none',legend.title=element_text(size=24,face='bold'),
              axis.text=element_text(size=20),axis.title=element_text(size=24),
              plot.title = element_text(size=36,hjust = 0.5))

  print(gp)
  ggsave(filename=paste0('./outputs/Trans_prot_compare/',x,'Transcriptome_Proteome_Compare.tiff'),plot=gp,device='tiff',height=10,width=16)

  #writeLines(high.prot,paste0('./outputs/gene_groups/',x,'Proteome_but_low_expression.txt'))
  #writeLines(high.trans,paste0('./outputs/gene_groups/',x,'High_expresssion_but_no_proteome.txt'))

  return(l)

}

samps=c('L2','L3','L4','L5_FG','L5_AMG','L5_MMG','L5_PMG','L5_HG')
l2=lapply(samps,omic.compare)



full=rbindlist(unlist(l2,recursive=F))
high.prot=full %>% filter(proteome=='Present') %>% filter(transcriptome_adjust<low.thresh) %>% .[['ha_geneid']]
high.trans=full %>% filter(proteome=='Absent') %>% filter(transcriptome_adjust>high.thresh) %>% .[['ha_geneid']]

writeLines(high.prot,'./outputs/gene_groups/Combined_Proteome_but_low_expression.txt')
writeLines(high.trans,'./outputs/gene_groups/Combined_High_expresssion_but_no_proteome.txt')


#### Make plot
gp=ggplot(data=full,aes(x=proteome,y=transcriptome_adjust))
gp=gp+geom_violin()
gp=gp+scale_fill_aaas()
gp=gp+labs(x='\nProteome Call',y='Expression Level (TPM)\n')
gp=gp+ggtitle('Transcriptomic and Proteomic Comparison')
gp=gp+scale_y_continuous(trans = log2_trans(),breaks=c(2,5,10,20,50,100,200,400,800,1500,3000))
gp=gp+geom_segment(aes(x=.5,xend=1.5,y=high.thresh,yend=high.thresh),linetype=2)
gp=gp+geom_segment(aes(x=1.5,xend=2.5,y=low.thresh,yend=low.thresh),linetype=2)

gp=gp+theme_bw()
gp=gp+theme(legend.position='none',legend.title=element_text(size=24,face='bold'),
            axis.text=element_text(size=20),axis.title=element_text(size=24),
            plot.title = element_text(size=36,hjust = 0.5))

print(gp)
ggsave(filename='./outputs/Trans_prot_compare/Combined_Transcriptome_Proteome_Compare.png',plot=gp,device='png',height=10,width=16)





# 
# 
# # ###################### Quantitative comparison #######################
# quant.omic.compare=function(x){
#  #x='L2'
#   l=list()
#  sub=merged.omics.sum %>% select(ha_geneid,matches(x),-matches('_C_'))
# 
# 
#  prot=colnames(sub)[grepl('prot',colnames(sub))]
#  trans=colnames(sub)[grepl('trans',colnames(sub))]
# 
#  sub=sub[get(prot)!=0]
# 
#  #sub$trans_adjust=sub[[trans]]+1
#  #sub$prot_adjust=sub[[trans]]+1
# 
#  sub2=sub %>% `colnames<-`(c('ha_geneid','transcriptome','proteome')) %>% mutate(sample=x)
#  l[[x]]=sub2
# 
#  ###extract outliers
#  #high.prot=sub %>% filter(get(prot)=='Present') %>% filter(get(trans)<low.thresh) %>% .[['ha_geneid']]
# 
#  #high.trans=sub %>% filter(get(prot)=='Absent') %>% filter(get(trans)>high.thresh) %>% .[['ha_geneid']]
# 
#  #### Make plot
#  gp=ggplot(data=sub,aes_string(x=prot,y=trans))
#  gp=gp+geom_point()
#  #gp=gp+scale_fill_aaas()
#  gp=gp+labs(x='\nProteome Call',y='Expression Level (TPM)\n')
#  gp=gp+ggtitle(paste0(x,' Transcriptomic and Proteomic Comparison'))
#  #gp=gp+scale_y_continuous(trans = log2_trans(),breaks=c(2,5,10,20,50,100,200,400,800,1500,3000))
#  gp=gp+scale_y_continuous(trans = pseudo_log_trans(base = 10),breaks=c(2,5,10,20,50,100,200,400,800,1500,3000))
#  gp=gp+scale_x_continuous(trans = pseudo_log_trans(base = 10),breaks=c(2,5,10,20,50,100,200,500,1000,1500,3000))
#  gp=gp+theme_bw()
#  gp=gp+theme(legend.position='none',legend.title=element_text(size=24,face='bold'),
#              axis.text=element_text(size=20),axis.title=element_text(size=24),
#              plot.title = element_text(size=36,hjust = 0.5))
# 
#  #print(gp)
#  ggsave(filename=paste0('./outputs/Trans_prot_compare/',x,'Transcriptome_Proteome_Compare_QUANT.tiff'),plot=gp,device='tiff',height=10,width=16)
# 
#  #writeLines(high.prot,paste0('./outputs/gene_groups/',x,'Proteome_but_low_expression.txt'))
#  #writeLines(high.trans,paste0('./outputs/gene_groups/',x,'High_expresssion_but_no_proteome.txt'))
# 
#  return(l)
# 
# }
# 
# 
# samps=c('L2','L3','L4','L5_FG','L5_AMG','L5_MMG','L5_PMG','L5_HG')
# l2=lapply(samps,quant.omic.compare)
# 
# 
# full=rbindlist(unlist(l2,recursive=F))
# full.gel=full[grepl('L5',sample)] %>% filter(proteome!=0)
# full.free=full[!grepl('L5',sample)] %>% filter(proteome!=0)
# 
# 
# # 
# ############## Histogram of values
# gp=ggplot(data=full.gel,aes(x=proteome))
# gp=gp+geom_histogram()
# print(gp)
# 
# 
# with(full.gel,cor.test(proteome,transcriptome))
# with(full.free,cor.test(proteome,transcriptome))
# 
# gp=ggplot(data=full.gel,aes(x=proteome,y=transcriptome))
# gp=gp+geom_point()
# #gp=gp+scale_fill_aaas()
# gp=gp+labs(x='\nProteome Call',y='Expression Level (TPM)\n')
# gp=gp+ggtitle('Global Quantitative L5 Comparison')
# gp=gp+geom_smooth(method='lm')
# gp=gp+scale_y_continuous(trans = pseudo_log_trans(base = 10),breaks=c(2,5,10,20,50,100,200,400,800,1500,3000))
# gp=gp+scale_x_continuous(trans = pseudo_log_trans(base = 10),breaks=c(.02,.04,.08,.2,.4,1,2,5,10,20,50,100,200,500,1000,1500,3000))#,limits = c(0,.2))
# gp=gp+theme_bw()
# gp=gp+theme(legend.position='none',legend.title=element_text(size=24,face='bold'),
#            axis.text=element_text(size=20),axis.title=element_text(size=24),
#            plot.title = element_text(size=36,hjust = 0.5))
# 
# print(gp)
# ggsave(filename=paste0('./outputs/Trans_prot_compare/Gel_Based_TransPort_Compare_Quant.png'),plot=gp,device='png',height=10,width=16)
# 
# 
# 
# gp=ggplot(data=full.free,aes(x=proteome,y=transcriptome))
# gp=gp+geom_point()
# #gp=gp+scale_fill_aaas()
# gp=gp+labs(x='\nProteome Call',y='Expression Level (TPM)\n')
# gp=gp+ggtitle('Global Quantitative L5 Comparison')
# gp=gp+geom_smooth(method='lm')
# gp=gp+scale_y_continuous(trans = pseudo_log_trans(base = 10),breaks=c(2,5,10,20,50,100,200,400,800,1500,3000))
# gp=gp+scale_x_continuous(trans = pseudo_log_trans(base = 10),breaks=c(2,5,10,20,50,100,200,500,1000,1500,3000))
# gp=gp+theme_bw()
# gp=gp+theme(legend.position='none',legend.title=element_text(size=24,face='bold'),
#             axis.text=element_text(size=20),axis.title=element_text(size=24),
#             plot.title = element_text(size=36,hjust = 0.5))
# 
# print(gp)
# ggsave(filename=paste0('./outputs/Trans_prot_compare/Transcriptome_Proteome_Compare_QUANT.tiff'),plot=gp,device='tiff',height=10,width=16)
# 
# 
# 
# 
