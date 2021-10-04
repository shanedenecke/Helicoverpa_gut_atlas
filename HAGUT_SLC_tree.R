#!/usr/bin/env Rscript  

####################################################  Import Libraries #################################################### 
shhh <- suppressPackageStartupMessages
shhh(library(dplyr)) ### maniuplating data
shhh(library(data.table)) ### alternative data.frame structure
shhh(library(ape)) ## useful for
shhh(library(ggtree)) ### Tree plotting package add on to ggplot2
shhh(library(ggplot2))### basic plotting package
shhh(library(ggstance)) ### horizontal plots for facet_plot
shhh(library(phytools)) ### midpoint root function


#################################################### Set working directory  #################################################### 
setwd('~/Dropbox/Shane/omics_projects/Ha_midgut_atlas')
#setwd('C:/Users/shane/Dropbox/Shane/omics_projects/Ha_midgut_atlas')

dir.create('./outputs/SLC_trees',showWarnings = F)

#################################################### Import Functions #################################################### 

### This function collapses nodes on the tree depending on what value you set
di2multi4node <- function (phy, tol = 0.5) {
  if (is.null(phy$edge.length)) {stop("the tree has no branch length")}
  if (is.na(as.numeric(phy$node.label[2]))){stop("node labels can't be converted to numeric values")}
  if (is.null(phy$node.label)){stop("the tree has no node labels")}
  ind <- which(phy$edge[, 2] > length(phy$tip.label))[as.numeric(phy$node.label[2:length(phy$node.label)]) < tol]
  n <- length(ind)
  if (!n) {return(phy)}
  foo <- function(ancestor, des2del) {
    wh <- which(phy$edge[, 1] == des2del)
    for (k in wh) {
      if (phy$edge[k, 2] %in% node2del) 
        foo(ancestor, phy$edge[k, 2])
      else phy$edge[k, 1] <<- ancestor
    }
  }
  node2del <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node2del) {next}
    foo(anc[i], node2del[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) {phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]}
  phy
}


### This funciton calcualtes standard error. Why is it not included in R
sterr=function(x) sd(x)/sqrt(length(x))

#### remove genes with below a certain level of average expression
low.remove=function(data,thresh=5,id='ha_geneid'){
  #data=trans.sum
  max_expression=split(data,data[[id]]) %>% lapply(.,function(x) max(x$average_expression))
  good_genes=max_expression[max_expression>thresh]
  return(data[get(id) %in% names(good_genes)])
}



#### label nodes with color bootstraps

treecol=function(tree){
  num=as.numeric(tree$node.label)
  cols=c()
  for(j in num){
    if(is.na(j)){cols=c(cols,'black')
    }else if(j>=80){cols=c(cols,'green1')
    }else if(j>=50){cols=c(cols,'cyan')
    }else{cols=c(cols,'red')}
  }
  return(cols)
}


####################################################  Import tree and collapse nodes #################################################### 
slc2.tree=read.tree('./inputs/transporters/SLC2_HelArm.nwk') %>% di2multi4node(.,as.numeric(50)) %>% midpoint.root()#root(.,outgroup='Helicoverpa_armigera___HelArm___110382520_XP_021198813.1___SLC_15_1')
slc22.tree=read.tree('./inputs/transporters/SLC22_HelArm.nwk') %>% di2multi4node(.,as.numeric(50)) %>% midpoint.root() #root(.,outgroup='Helicoverpa_armigera___HelArm___110377832_XP_021192542.1___SLC_37_1')

slc2.tree$tip.label=gsub("Helicoverpa_armigera___HelArm___([0-9]+)_.+$",'LOC\\1',slc2.tree$tip.label) ### subset out Geneids from long SLC names
slc22.tree$tip.label=gsub("HelArm___([0-9]+)_.+$",'LOC\\1',slc22.tree$tip.label) ### subset out Geneids from long SLC names

slc2.table=as_tibble(slc2.tree) %>% data.table() ### get table of phylogeny with nodes tip labels etc.
slc22.table=as_tibble(slc22.tree) %>% data.table() ### get table of phylogeny with nodes tip labels etc.

#rt=table[grepl('1167',label)]$label ## get root


#################################################### Process expression values for SLC_table #################################################### 
trans.raw=fread('./inputs/transporters/HelArm_SLC_Transporters_Clean.csv') %>% 
  select(-ha_proteinid) %>% filter(fam %in% c('SLC_2','SLC_22','SLC_15','SLC_37')) %>% 
  mutate(ha_geneid=paste0('LOC',ha_geneid))
tpm.melted=fread('./outputs/omics_clean/Ha_Transcripomics_TPM_Replicate_Melt.csv')

trans.full=tpm.melted %>% merge(trans.raw,by='ha_geneid') %>% 
  filter(Tissue %in% c("L2_MG","L3_MG","L4_MG")) %>% 
  rename(label=ha_geneid) %>% 
  mutate(TPM=TPM+1)

trans.sum=trans.full %>% group_by(label,Tissue,fam) %>% 
  summarize(average_expression=mean(TPM),
            low=mean(TPM)-sterr(TPM),
            high=mean(TPM)+sterr(TPM)) %>% data.table %>% low.remove(.,2,'label')




#################################################### Subset out each SLC family and then prune trees #################################################### 
slc2.sum=trans.sum[fam=='SLC_2']
slc2.tree=drop.tip(slc2.tree,slc2.tree$tip.label[!(slc2.tree$tip.label %in% slc2.sum$label)])


slc22.sum=trans.sum[fam=='SLC_22']
slc22.tree=drop.tip(slc22.tree,slc22.tree$tip.label[!(slc22.tree$tip.label %in% slc22.sum$label)])



#################################################### Generate Heatmaps #################################################### 
meanNorm=function(data){
  x=data$TPM
  minus=(x-mean(x))
  #trans=minus/(sd(x))
  trans=minus/(max(x)-min(x))
  data$normTPM=trans
  return(data)
}
heat=select(trans.full,fam,label,Tissue,TPM) %>% filter(fam=='SLC_2') %>%
  group_by(fam,Tissue,label) %>% summarize(TPM=mean(TPM)) %>% filter(label %in% slc2.tree$tip.label) %>%
  mutate(TPM=log2(TPM))

heat2=split(heat,heat$label) %>% lapply(meanNorm) %>% rbindlist()


# gp=ggplot(data=heat2,aes(x=Tissue,y=label,fill=normTPM))
# gp=gp+scale_fill_gradient(low='blue',high='red')
# gp=gp+geom_tile(color='black')
# gp=gp+labs(x='\nSpecies 1',y='Species 2\n')
# gp=gp+theme_bw()
# gp=gp+theme(legend.title=element_text(size=18,face='bold'),
#             legend.text = element_text(size=12,face='bold'),
#             axis.text = element_text(size=10),
#             axis.text.x=element_text(angle=90),
#             axis.title=element_text(size=24),
#             plot.title = element_text(size=36,hjust = 0.5))
# print(gp)

heat3=heat2 %>% pivot_wider(id_cols=label,names_from=Tissue,values_from = normTPM) %>% data.frame()
colnames(heat3)=gsub('_',' ',colnames(heat3))
rownames(heat3)=heat3$label
heat3$label=NULL
heat3=as.matrix(heat3)




#################################################### Get GFF file to annotate scaffold #################################################### 
gtf=fread('../omics_ref/Helicoverpa_armigera/ncbi/harmi_ABC123_mod.gtf')
gtf$V10=gsub('.+(LOC[0-9]+).+','\\1',gtf$V9)
gtf.reduce=distinct(gtf,V1,V10) %>% filter(V10 %in% unique(slc2.tree$tip.label))
scaff.col=c()
for(i in unique(slc2.tree$tip.label)){
  scaff.col=c(scaff.col,ifelse(gtf.reduce[V10==i]$V1=='NW_018395529.1','orange4','black'))
  }

#################################################### SLC2 Tree #################################################### 
gp=ggtree(slc2.tree,size=1)
gp=gp+geom_tippoint(size=4,color=scaff.col)
gp=gp+geom_nodepoint(size=3,col=treecol(slc2.tree))
gp=gp+xlim_tree(4)
gp=gp+labs(x='Substitutions / Site')
gp=gp+geom_tiplab(fontface='bold',hjust=-.3,align=T,color=scaff.col)
gp=gp+theme(legend.text=element_text(size=30,face='bold'),legend.title = element_text(size=36,face='bold'))
gp=gp+theme_tree2()


#  gp2=facet_plot(gp,data=slc2.sum,panel = "Midgut Expression", geom=geom_barh,
#                 mapping=aes(x=average_expression,fill=Tissue),
#                 width=.8,position=position_dodgev(),stat='identity',color=NA)
#  
#  gp3=facet_plot(gp2,panel='Midgut Expression',geom=geom_hline,data=slc2.sum,
#                mapping=aes(yintercept = y+.5),linetype=2)
# 
# ### Ignore warning about aestetic
#  gp4=facet_plot(gp3,data=slc2.sum,panel = "Midgut Expression", geom=geom_errorbarh,
#                 aes(x=average_expression,xmin=low,xmax=high,group=interaction(label,Tissue)),
#                 position=position_dodgev(),width=.8,color='black')
# 



gp4=gheatmap(gp, heat3, offset=.8, width=0.6, 
         colnames=T, legend_title="Gene Expression \n (Mean Normalized)",
         low='blue',high='red',color='black',colnames_position = 'top',font.size = 6)

gp4=gp4+labs(x='Total Expression of gene (log2(TPM))')
#gp4=gp4+theme_tree2()
gp4=gp4+labs(x='Substitutions / Site')
gp4=gp4+theme(strip.text = element_text(size=24),
              legend.text=element_text(size=14),legend.title=element_text(size=20,face='bold'),
              axis.text=element_text(size=14),axis.title=element_text(size=14,hjust=.8))

print(gp4)



ggsave(filename='./outputs/SLC_trees/SLC2_Life_Stage.pdf',plot=gp4,height=9,width=12,device='pdf')









#################################################### SLC22 Tree #################################################### 
gp=ggtree(slc22.tree,size=1)
gp=gp+geom_tippoint(size=4)
gp=gp+xlim_tree(5)
gp=gp+geom_tiplab(fontface='bold',hjust=-.3,align=T)
gp=gp+theme(legend.text=element_text(size=30,face='bold'),legend.title = element_text(size=36,face='bold'))

gp2=facet_plot(gp,data=slc22.sum,panel = "Midgut Expression", geom=geom_barh,
               mapping=aes(x=average_expression,fill=Tissue),
               width=.8,position=position_dodgev(),stat='identity',color=NA)

gp3=facet_plot(gp2,panel='Midgut Expression',geom=geom_hline,data=slc22.sum,
               mapping=aes(yintercept = y+.5),linetype=2)

#### Ignore warning about aestetic
gp4=facet_plot(gp3,data=slc22.sum,panel = "Midgut Expression", geom=geom_errorbarh,
               aes(x=average_expression,xmin=low,xmax=high,group=interaction(label,Tissue)),
               position=position_dodgev(),width=.8,color='black')


gp4=gp4+labs(x='Total Expression of gene (log2(TPM))')

gp4=gp4+theme_tree2()
gp4=gp4+theme(strip.text = element_text(size=24),
              legend.text=element_text(size=14),legend.title=element_text(size=20,face='bold'),
              axis.text=element_text(size=14),axis.title=element_text(size=14,hjust=.8))

gp4

ggsave(filename='./outputs/SLC_trees/SLC22_Life_Stage.pdf',gp4,height=9,width=12,device='pdf')

