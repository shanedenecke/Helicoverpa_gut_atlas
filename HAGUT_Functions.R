reduce.genes=function(x,column,annot){
  l=list()
  for(i in names(x)){
    el=x[[i]]
    l[[i]]=el[[column]]
  }
  return(annot[get(column) %in% Reduce(intersect,l)])
}


freq.table=function(x){
  total=x$Ha_gene %>% unique() %>% length()
  l=list()
  for(i in unique(x$code)){
    sub=x[code==i]
    num.genes=unique(sub$Ha_gene) %>% length()
    freq=num.genes/total
    l[[i]]=data.table(code_name=i,code_count=num.genes,frequency=freq,total=total)
  }
  return(rbindlist(l))
} 

vlepo.enrich=function(trial,omic){
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
  enriched= combined %>% arrange(fdr) %>% filter(fdr<1e-02) %>% filter(test_N>10) %>% select(-frequency,-code_count) %>% 
    mutate(test_freq=test_N/test_total, ref_freq=ref_N/ref_total) %>% 
    select(code_name,fischP,fdr,everything()) %>% data.table()
  return(enriched)
}

annotate_enrichment=function(enrich.table,key){
  l=list()
  for(i in names(enrich.table)){
    enrich=enrich.table[[i]]
    annot=merge(enrich,key,by='code_name',all.y=F,sort=F)
    l[[i]]=(annot[!duplicated(annot)])
  } 
  return(l)
}

gene.subset=function(geneset,database){
  return(database[Ha_gene %in% geneset])
}




gene.extract=function(x,column){
  return(x[[column]])
}
