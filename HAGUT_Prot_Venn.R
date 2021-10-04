#!/usr/bin/env Rscript

######################################## Import Packages ########################################

shhh=suppressPackageStartupMessages
shhh(library(data.table))
shhh(library(tidyverse))
shhh(library(VennDiagram))
shhh(library(UpSetR))

######################################## Directories ########################################
setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
source('./scripts/HAGUT_Functions.R')
dir.create('./outputs/Venn_diagrams',showWarnings = F) 


######################################## Import Data ########################################

prot=fread('./outputs/proteomics_clean/Final_Proteomics_Presence.csv') %>% 
  #select(ha_geneid,matches('prot')) %>% select(-contains('std')) %>% 
  filter_at(vars(matches('prot')),any_vars(.=='Present'))


######################################## Define Functions  ########################################
prot.unique=function(x,data){
  uni=data %>% filter_at(vars(contains(x)),all_vars(.=='Present')) %>% 
    filter_at(vars(-contains(x), -contains('ha_geneid')), all_vars(.=='Absent')) %>%
    data.table()
  return(uni$ha_geneid)
}

prot.common=function(data){
  uni=data %>% filter_at(vars(contains('prot')),all_vars(.=='Present'))
  return(list(uni$ha_geneid))
}


######################################## Split Data into Lists ########################################
samps=colnames(select(prot,-ha_geneid))

l=lapply(samps,function(x) prot[get(x)=='Present']$ha_geneid)
names(l)=samps
names(l)=gsub('_prot','',names(l))

spatial=l[(names(l) %in% c("L5_FG",'L5_AMG','L5_MMG','L5_PMG','L5_HG'))]
life.stage=l[names(l) %in% c("L3",'L2','L4')]

######################################## Make Venn Diagrams ########################################

###spatial venn
spatial.ven=venn.diagram(spatial,filename=NULL,main='Spatial Comparison',cex=1.5,fill=c('coral','cadetblue','gold','lightgreen','grey50'))
ggsave(spatial.ven,file='./outputs/Venn_diagrams/Spatial_Comparison.tiff',device='tiff',width=15,height=10)

pdf('./outputs/Venn_diagrams/Upside_Plot_spatial.pdf',width=20,height=12)
upset.spatial=upset(fromList(spatial),keep.order=T,order.by = "freq",mainbar.y.label = "Number of Common Proteins", 
                    sets=c("L5_HG", "L5_PMG", "L5_MMG", "L5_AMG", "L5_FG"),
                    sets.x.label = "Total Proteins Detected in Each Sample", 
                    text.scale = c(4, 2, 2, 2, 2,2),
                    line.size=2,
                    point.size=8,
                    sets.bar.color='black',
                    main.bar.color='black',
                    #group.by='sets',
                    )
print(upset.spatial)
dev.off()


###Midgut only compartments
midgut.only=spatial[grepl('MG',names(spatial))]
midgut.only.ven=venn.diagram(midgut.only,filename=NULL,main='Midgut Spatial Comprison',cex=1.5,fill=c('coral','cadetblue','gold'))
ggsave(midgut.only.ven,file='./outputs/Venn_diagrams/Only_midgut.tiff',device='tiff',width=15,height=10)

###Life stage comprison
life.ven=venn.diagram(life.stage,filename=NULL,main='Life_stage Comparison',cex=1.5,fill=c('coral','cadetblue','gold'))
ggsave(life.ven,file='./outputs/Venn_diagrams/Life_stage.tiff',device='tiff',width=10,height=10)



######################################## Extract gene lists ######################################## 

### subset spatial and stage data tables
spatial.sub=prot %>% select(ha_geneid,L5_FG_prot,L5_AMG_prot,L5_MMG_prot,L5_PMG_prot,L5_HG_prot)
stage.sub=prot %>% select(ha_geneid,L2_prot,L3_prot,L4_prot)
midgut.sub=prot %>% select(ha_geneid,L5_AMG_prot,L5_MMG_prot,L5_PMG_prot)


##extract unique and common protein subset functions
spatial.list.unique=lapply(c('FG','AMG','MMG','PMG','HG'),prot.unique,spatial.sub)
stage.list.unique=lapply(c('L2','L3','L4'),prot.unique,stage.sub)
spatial.list.common=prot.common(spatial.sub)
stage.list.common=prot.common(stage.sub)
midgut.common=list(spatial.sub[L5_AMG_prot=='Present' & L5_MMG_prot=='Present' & L5_PMG_prot=='Present' & L5_FG_prot=='Absent' & L5_HG_prot=='Absent']$ha_geneid)


######################################## Combine and Write final proteome lists ######################################## 
final.list=c(spatial.list.unique,stage.list.unique,spatial.list.common,stage.list.common,midgut.common)
names(final.list)=c('FG_unique','AMG_unique','MMG_unique','PMG_unique','HG_unique',
                    'L2_unique','L3_unique','L4_unique',
                    'Compartment_common','Stage_common','Midgut_common')


for(i in names(final.list)){
  writeLines(final.list[[i]],con=paste0('./outputs/gene_groups/',i,'_proteome.txt'))
}


for(i in list.files('./',full.names = T)){if(grepl('log',i)){file.remove(i)}}


