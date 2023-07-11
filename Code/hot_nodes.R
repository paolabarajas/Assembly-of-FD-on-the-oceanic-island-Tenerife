hot.nodes <- function(Sample, Tree, runs = 999, signif = 1.96, min.size = 10, positive = TRUE, verbose = TRUE){
  
  if(names(table(names(Sample)==Tree$tip.label))[[1]]=="TRUE") {
    1
  } else
  {stop("names in Sample are not in the same order as in the phylogeny")}
  
  Table <- data.frame(Tree$tip.label,1:length(Tree$tip.label))
  Tips_frac <-vector(mode="list",length=length(Tree$tip.label)-1)
  length(Tree$tip.label)+1->Start
  (length(Tree$tip.label)-1)+(length(Tree$tip.label))->End
  
  # To obtain species subtended in node i
  
  for(i in Start:End) {
    getDescendants(Tree, i, curr=NULL)->tips.i
    Table.i <- Table[Table[,2]%in%tips.i,]
    Table.i[,1]->Tips_frac[[i-(Start-1)]]
  }
  
  # To extract subtrees
  
  Trees_store <- vector(mode="list", length=length(Tips_frac))
  for(j in 1:length(Tips_frac)) {
    Trees_store[[j]]<-drop.tip(Tree,Tree$tip.label[Tree$tip.label%in%Tips_frac[[j]]==FALSE])
  }
  
  L_linajes <- vector(mode="numeric",length=length(Trees_store))
  for(z in 1:length(Trees_store)){
    Trees_store[[z]]->Tree.i
    length(Tree.i$tip.label)->L_linajes[z]
  }
  
  ####
  
  if(positive == TRUE) {
    
    Sample_1 <- Sample[Sample>0]
    
  }
  
  if(positive == FALSE) {
    
    Sample_1 <- Sample[Sample==0]
    Sample_1[]<-1
    Sample[Sample==1]<-2
    Sample[Sample==0]<-1
    Sample[Sample==2]<-0
    
  }
  
  Tree$tip.label[Tree$tip.label%in%names(Sample_1)]->Tree_Sample_1
  Tree$tip.label[Tree$tip.label%in%names(Sample_1)==FALSE]->Tree_Sample_0
  
  Sample_All<-c(rep(1,length(Tree_Sample_1)),rep(0,length(Tree_Sample_0)))
  names(Sample_All)<-c(Tree_Sample_1,Tree_Sample_0)
  
  SES_Sample <- rep(NA,length(Trees_store))
  Long_Lin <- rep(NA,length(Trees_store))
  Usadas <- rep(NA,length(Trees_store))
  
  for(i in 2:length(Trees_store)){
    
    Trees_store[[i]]->Tree.i
    length(Tree.i$tip.label[Tree.i$tip.label%in%Tree_Sample_1])->En_uso
    En_uso->Usadas[i]
    Store_NULL <- rep(NA,runs)
    
    for(j in 1:runs){
      
      sample(Sample_All)->Sample_NULL
      names(Sample_NULL)<-names(Sample_All)
      Sample_NULL_1<-Sample_NULL[Sample_NULL>0]
      
      length(Tree.i$tip.label[Tree.i$tip.label%in%names(Sample_NULL_1)]) -> Store_NULL[j]
      
    }
    
    length(Tree.i$tip.label)->Long_Lin[i]
    
    (En_uso-mean(Store_NULL))/sd(Store_NULL)->SES_Sample[i]
    
    if(verbose == TRUE){
      print(paste("node ",i-1," out of ",(length(Trees_store)-1)," completed",sep=""))} else
      {1}
  }
  
  data.frame(SES_Sample,Long_Lin,Usadas,1:length(Trees_store))->SES_Sample
  
  if(positive == TRUE) {
    colnames(SES_Sample)<-c("SES","Size","Usable","Node")
  }
  if(positive == FALSE) {
    colnames(SES_Sample)<-c("SES","Size","Non-usable","Node")
  }
  
  Tomados_all <- SES_Sample[-1,]
  
  Tomados_all <- Tomados_all[,c(4,2,3,1)]
  
  SES_Sample[SES_Sample[,1] > signif & SES_Sample[,2] > (min.size-1) ,]->Tomados
  
  Tomados <- Tomados[-1,]
  
  Tomados <- Tomados[,c(4,2,3,1)]
  
  Result<-list(Tomados_all, Tomados)
  
  if(positive == TRUE) {
    names(Result)<-c("All nodes","Hot nodes")
  }
  if(positive == FALSE) {
    names(Result)<-c("All nodes","Cold nodes")
  }
  return(Result)
}