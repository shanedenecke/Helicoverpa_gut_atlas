#!/usr/bin/env Rscript

######################################## import packages ######################################## 
shhh <- suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(scales))
shhh(library(ggsci))

###Set working directory
setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
dir.create('./outputs/pH',showWarnings = F)

######################################################################Import genes and summarized transcriptomes###################################


tpm.bySample=fread('./outputs/omics_clean/Ha_proteomics_transcriptomics_Replicates.csv')
ph=fread('./inputs/acid_base/HelArm_pH_genes.csv') %>% select(-ha_proteinid) %>% mutate(ha_geneid=paste0('LOC',ha_geneid))
tpm.melted=fread('./outputs/omics_clean/Ha_Transcripomics_TPM_Replicate_Melt.csv')
tpm.sum=fread('./outputs/omics_clean/Ha_summarized_proteomics_transcriptomics.csv')


######################################################################Functions ###################################
geom.series=function(start,stop){
  final=start
  if(start==0){
    i=10; final=c(final,i)
    }else{
      i=start
    }
  while(i<stop){
    i=i*2
    final=c(final,i)
  }
  last=final[length(final)]*2
  #final=c(final,last)
  final=sapply(final,signif,digits=2)
  return(final)
}

shane.transpose=function(dt,newcol){
  new.rows=colnames(dt)[2:length(colnames(dt))]
  a=transpose(dt)
  colnames(a)=as.character(a[1,])
  a=a[-1,]
  b=mutate(a,newcol=new.rows) %>% select(newcol,everything()) %>% data.table()
  return(b)
}

zs=function(x){
  m=mean(x$TPM)
  dev=sd(x$TPM)
  x$z_score=(x$TPM-m)/dev
  return(x)
}



######################################################################merge pH annotations with omics data###################################
pH.bySample=tpm.melted %>% merge(ph,by='ha_geneid') %>% 
  filter(Diet_type!='Artificial') %>% mutate(Tissue=factor(Tissue,levels=c('L5_FG','L5_AMG','L5_MMG','L5_PMG','L5_HG')))
fwrite(pH.bySample,'./outputs/pH/pH_gene_expression_bySample.csv')
  
ratio.sum=pH.bySample %>% group_by(ha_geneid) %>% 
  summarise(AMG_ratio=mean(TPM[Tissue=='L5_AMG'])/mean(TPM[Tissue=='L5_PMG']),
            total_exp=sum(TPM),
            Category=unique(Category),
            Name=unique(Name)) %>%
  data.table()

######################################################################Combined vATPase CAH Ratio###################################

all.sum=ratio.sum[Category=='CAH' | Category=='vATPase' | Category=='SLC_9' | Category=='SLC_26'] %>% filter(total_exp>25) 
gp=ggplot(all.sum,aes(x=Category,y=AMG_ratio,shape=Category))
gp=gp+ stat_summary(geom = "errorbar", colour = "black", size = 1, fun.data = mean_se,  width = 0.2)
gp=gp+stat_summary(fun = "mean", size=1, geom = "crossbar", colour = "black",width = 0.3)
gp=gp+geom_point(size=7,position=position_jitterdodge(jitter.width=.3,seed=2),aes(shape=Category,fill=Category),color='black',alpha=.6)
#gp=gp+geom_dotplot(binaxis = "y", stackdir = "center", binwidth = .2, colour = "black",dotsize = 1,fill='grey50')#,aes(fill=Name),position=position_jitter(width=.03,height=.03,seed=4))
gp=gp+scale_fill_aaas()
gp=gp+scale_shape_manual(values=c(21,22,23,24))
gp=gp+labs(x='\nGene Group',y='AMG Expression / PMG Expression')
#gp=gp+lims(y=c(0,5.1))
gp=gp+geom_hline(yintercept = 1,linetype=2)
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
print(gp)
ggsave('./outputs/pH/CAH_vATPase_ratio.png',gp,device='png',height=8,width=8)




######################################################################carbonic anhydrase full######################################################################
cah.all=pH.bySample[Category=='CAH'] %>% filter(Name %in% all.sum$Name) %>% mutate(TPM=TPM+1)


gp=ggplot(cah.all,aes(x=Tissue,y=TPM,group=Name,color=Name))
#gp=gp+facet_grid(.~Category,scales='free_y')
gp=gp+stat_summary(geom='point',fun='mean',size=3)
gp=gp+stat_summary(geom='line',fun='mean',size=1)
gp=gp+stat_summary(fun.data = mean_se, geom = "errorbar",width=.1)
gp=gp+scale_color_futurama()
gp=gp+ scale_y_continuous(trans = log2_trans())
gp=gp+ggtitle('Carbonic Anhydrase expression along the midgut')
gp=gp+labs(x='Diet Condition',y='Transcripts per Million')
gp=gp+geom_hline(yintercept = 1,linetype=2)
gp=gp+theme_bw()
gp=gp+theme(axis.title=element_text(size=20,face='bold'),
            panel.grid=element_blank(),
            axis.text.x=element_text(size=14,face='bold',hjust=1,angle=45),
            axis.text.y=element_text(size=14,face='bold'),
            rect=element_rect(fill=NULL),
            #legend.position='none',
            strip.text = element_text(size=20,face='bold'),
            plot.title = element_text(size=28,face='bold',hjust = 0.5))
print(gp)   
ggsave('./outputs/pH/Carbonic_anhydrase_spatial.png',gp,device='png',height=8,width=8)







######################################################################vATPase full######################################################################
va.all=pH.bySample[Category=='vATPase'] %>% filter(Name %in% all.sum$Name) %>% filter(Name!='Vha100-2') %>% mutate(TPM=TPM+1)

gp=ggplot(va.all,aes(x=Tissue,y=TPM,group=Name),color='black',alpha=.6)
gp=gp+stat_summary(geom='point',fun='mean',size=3)
gp=gp+stat_summary(geom='line',fun='mean',size=1)
gp=gp+stat_summary(fun.data = mean_se, geom = "errorbar",width=.1)
#gp=gp+scale_color_grey()
#gp=gp+lims(y=c(20,20000))
gp=gp+ scale_y_continuous(trans = log2_trans(),breaks = geom.series(20,20000))
gp=gp+ggtitle('vATPase expression along the midgut')
gp=gp+labs(x='Diet Condition',y='Transcripts per Million')
gp=gp+theme_bw()
gp=gp+theme(axis.title=element_text(size=20,face='bold'),
            panel.grid=element_blank(),
            axis.text.x=element_text(size=14,face='bold'),
            axis.text.y=element_text(size=14,face='bold'),
            rect=element_rect(fill=NULL),
            #legend.position='none',
            strip.text = element_text(size=20,face='bold'),
            plot.title = element_text(size=28,face='bold',hjust = 0.5))
print(gp)   

ggsave('./outputs/pH/vATPase_spatial.png',gp,device='png',height=8,width=8)


################################################ Other transporters######################################################################
slc.sum=ratio.sum %>% filter(grepl('SLC',Category)) %>% filter(!Category=='SLC_2') %>% filter(total_exp>20)

gp=ggplot(slc.sum,aes(x=Category,y=AMG_ratio))
gp=gp+ stat_summary(geom = "errorbar", colour = "black", 
                    fun.data = mean_se, width = 0.1,size=2)
gp=gp+stat_summary(fun = "mean", geom = "crossbar", colour = "black", 
                   width = 0.2,size=1)
gp=gp+geom_dotplot(binaxis = "y", 
                   stackdir = "center", binwidth = .5, colour = "black", 
                   position = position_dodge(width = 0.5), dotsize = 1,fill='grey50')
gp=gp+ggtitle('Ion transporter expression Enrichement')
gp=gp+labs(x='Diet Condition',y='Anterior Midgut Enrichment')
gp=gp+geom_hline(yintercept = 1,linetype=2)
gp=gp+theme_bw()
gp=gp+theme(axis.title=element_text(size=20,face='bold'),
            panel.grid=element_blank(),
            axis.text.x=element_text(size=14,face='bold',hjust=1,angle=45),
            axis.text.y=element_text(size=14,face='bold'),
            rect=element_rect(fill=NULL),
            legend.position='none',
            strip.text = element_text(size=20,face='bold'),
            plot.title = element_text(size=28,face='bold',hjust = 0.5))
print(gp)
ggsave('./outputs/pH/SLC_ion_ratio.png',gp,device='png',height=8,width=8)




slc=pH.bySample[grepl('SLC',Category)] %>% filter(!(Category=='SLC_2')) %>% mutate(TPM=1+TPM)

gp=ggplot(slc,aes(x=Tissue,y=TPM,group=Name,color=Category,shape=Category))
gp=gp+facet_grid(.~Category,scales = 'free_y')
gp=gp+stat_summary(geom='point',fun='mean',size=3)
gp=gp+stat_summary(geom='line',fun='mean',size=1)
gp=gp+stat_summary(fun.data = mean_se, geom = "errorbar",width=.1)
gp=gp+scale_color_aaas()
gp=gp+ scale_y_continuous(trans = log2_trans(),breaks = geom.series(0,max(slc$TPM)))
gp=gp+ggtitle('SLC expression along the midgut')
gp=gp+labs(x='Diet Condition',y='Transcripts per Million')
gp=gp+geom_hline(yintercept = 1,linetype=2)
gp=gp+theme_bw()
gp=gp+theme(axis.title=element_text(size=20,face='bold'),
            panel.grid=element_blank(),
            axis.text.x=element_text(size=14,face='bold',hjust=1,angle=45),
            axis.text.y=element_text(size=14,face='bold'),
            rect=element_rect(fill=NULL),
            legend.position='none',
            strip.text = element_text(size=20,face='bold'),
            plot.title = element_text(size=28,face='bold',hjust = 0.5))
print(gp)   
ggsave('./outputs/pH/Ion_transporters_total.png',width=18,height=9)



################################################vATPase SLC9######################################################################
va.slc9=pH.bySample[Category=='vATPase' | Category=='SLC_9'] %>% filter(Name %in% va.all$Name | Name %in% c('Nhe1','NHA1a'))
va.slc9=rbindlist(lapply(split(x=va.slc9,by=c('Replicate','Name')),zs))
va.slc9.sum=group_by(va.slc9,Name,Tissue) %>% summarize(m=mean(z_score)) %>% data.table() 

gp=ggplot(va.slc9,aes(x=Tissue,y=z_score,group=Name,color=Category))
gp=gp+stat_summary(geom='point',fun='mean',size=3,alpha=.3)
gp=gp+stat_summary(geom='line',fun='mean',size=2,alpha=.3)
gp=gp+stat_summary(fun.data = mean_se, geom = "errorbar",width=.1,alpha=.3)
gp=gp+scale_color_manual(values=c('blue','red'))
#gp=gp+ scale_y_continuous(trans = log2_trans())
#gp=gp+ggtitle('SLC expression along the midgut')
gp=gp+labs(x='\nCompartment',y='Normalized Z score (TPM)\n')
gp=gp+geom_hline(yintercept = 0,linetype=2,alpha=.3)
gp=gp+theme_bw()
gp=gp+theme(axis.title=element_text(size=20,face='bold'),
            panel.grid=element_blank(),
            axis.text.x=element_text(size=20,face='bold'),
            axis.text.y=element_text(size=20,face='bold'),
            legend.text=element_text(size=14),legend.title=element_text(size=20,face='bold'),
            rect=element_rect(fill=NULL),
            strip.text = element_text(size=20,face='bold'),
            plot.title = element_text(size=28,face='bold',hjust = 0.5))
print(gp)   
ggsave('./outputs/pH/vATPase_SLC9.png',width=12,height=6)




