#' return a list by BC
#' @param mdir
#' @param pathwaynames
#' @export 
BC<-function(mdir=NULL,pathwaynames=NULL){
	library(igraph)
	library(KEGGgraph)
	betweennesslist<-NULL
	pathwayID<-NULL
	for(i in 1:length(pathwaynames)){
		mapkpathway<-try(parseKGML(paste(mdir,pathwaynames[[i]],sep="/")),TRUE)
		mapkG<- KEGGpathway2Graph(mapkpathway, expandGenes=F)
		g<-igraph.from.graphNEL(mapkG)
		bet<-betweenness(g)
		nodlist<-NULL
		nodeslist<-NULL
		nod<-nodes(mapkpathway)
		for(j in 1:length(bet)){
			nodname<-names(bet[j])
			genename<-nod[[nodname]]@name
    
			for(jj in 1:length(genename)){
				betweenness<-rep(bet[nodname],length(genename))
			}
			nodlist<-t(rbind(genename,betweenness))
			nodeslist<-rbind(nodeslist,nodlist)
		}
		betness<-nodeslist[,2]
		betness<-1+as.numeric(betness)
		names(betness)<-nodeslist[,1]
		name<-names(betness)
		name<-strsplit(as.character(name),"hsa:")
		name<-do.call(rbind,name)
		names(betness)<-name[,2]
		betweennesslist<-c(betweennesslist,list(betness))
		pathwayID<-c(pathwayID,mapkpathway@pathwayInfo@name)
	}	
	names(betweennesslist)<-pathwayID
	return(betweennesslist)
}
########################################################

#' return a list by SP
#' @param mdir
#' @param pathwaynames
#' @export 
SP<-function(mdir=NULL,pathwaynames=NULL){
	library(igraph)
	library(KEGGgraph)
	nodeslist<-NULL
	pathwayID<-NULL
	for(i in 1:length(pathwaynames)){
		mapkpathway<-try(parseKGML(paste(mdir,pathwaynames[[i]],sep="/")),TRUE)
		mapkG<- KEGGpathway2Graph(mapkpathway, expandGenes=T)
		nodes<-nodes(mapkG)
		nodeslist<-c(nodeslist,nodes)
	}
	specific_number<-table(unlist(nodeslist))
	wfg<-specific_number
	wfg<-as.data.frame(wfg)
	nodenames<-wfg$Var1
	nodenames<-strsplit(as.character(nodenames),"hsa:")
	nodenames<-do.call(rbind,nodenames)
	wfg<-wfg[,2]
	names(wfg)<-nodenames[,2]
	return(wfg)
}

##############################################################

#' return a list by IF
#' @param mdir
#' @param pathwaynames
#' @param DE
#' @param ALL
#' @export 
IF<-function(mdir=NULL,pathwaynames=NULL,DE=NULL,ALL=NULL){
	library(igraph)
	library(KEGGgraph)
	all<-paste('hsa:',ALL,sep="")
	DE<-names(DE)
	DE<-paste('hsa:',DE,sep="")
	edge<-NULL
	nodeslist<-NULL
	for(i in 1:length(pathwaynames)){
		mapkpathway<-try(parseKGML(paste(mdir,pathwaynames[[i]],sep="/")),TRUE)
		mapkG<-KEGGpathway2Graph(mapkpathway,expandGenes=TRUE)
		node<-nodes(mapkG)
		edL<-edgeData(mapkG)
		circsp<-strsplit(as.character(names(edL)),"\\|")
		geneName<-do.call(rbind,circsp)
		edge<-rbind(edge,geneName)
		nodeslist<-c(nodeslist,node)
	}
	nodeslist<-nodeslist[!duplicated(nodeslist)]
	e<-unique.matrix(edge)
	g<-graph_from_edgelist(e)
	mapk<-igraph.to.graphNEL(g)
	nodes<-nodes(mapk)
	neighborhoods<-ego(g,1,nodes,"out")
	inter<-function(X){
		nde<-length(intersect(names(X),DE))
	}
	nDE<-lapply(neighborhoods,inter)
	nDE<-as.numeric(as.matrix(nDE))
	signodes<-setdiff(nodeslist,nodes)
	nDEsn<-rep(0,length(signodes))
	nDE<-c(nDE,nDEsn)
	names(nDE)<-c(nodes,signodes)
	nodenames<-strsplit(as.character(c(nodes,signodes)),"hsa:")
	nodenames<-do.call(rbind,nodenames)
	names(nDE)<-nodenames[,2]
	return(nDE)
}

##############################################################

#' return a list by IF
#' @param wfg
#' @param nDE
#' @export 
wi<-function(wfg=NULL,nDE=NULL){
	DEdegree<-nDE[names(wfg)]
	wig<-DEdegree/wfg
	write.csv(wig,"./GWSPIA/w_8671.csv")
	wig<-1+((wig-min(wig))/(max(wig)-min(wig)))
	return(wig)
}





