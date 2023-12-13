##' We built multilayer networks and calculated centralities based on the
##' approach of:
##' 
##'@references
##'De Domenico, M., Porter, M. A., & Arenas, A. (2015). MuxViz: a tool for 
##'multilayer analysis and visualization of networks. Journal of Complex 
##'Networks, 3(2), 159-176.
##'  \href{https://doi.org/10.1093/comnet/cnu038}{}
##'@export
##'
##'
##'For multilayer networks, pairwise connections between nodes in monolayer 
##'networks will now change based on sequential thresholds of 10 days
##'of LHC (inundation frequencies in days per year).
##'
##' @param yearly_graphLists list object containing sequence of SC and FC graph objects
##' @param Layers numeric vector describing the number of layers 
##' @param Nodes numeric vector describing the number of replica nodes (same 
##' number of nodes in each layer)
##' @param OmegaParameter (layerCouplingStrength) numeric vector describing the 
##' interlayer link weights. The default in muxviz is 1 (unweighted). 
##' @param MultisliceType (networkOfLayersType) numeric vector describing the 
##' multilayer type: categorical (bidirectional links between layers);
##' ordinal (unidirectional links between layers. For temporal networks, 
##' following the order t1->t2->tn)


#3. Built multilayer networks and calculate versatilities
get_multicentralities<- function(yearly_graphLists,
                                 L,Nodes,
                                 layerCouplingStrength,
                                 networkOfLayersType){
  
  
  nodeTensor_SC<- list()
  nodeTensor_FC<- list()
  
  #1. Built network tensor 
  for(l in 1:L){
    nodeTensor_SC[[l]] <- igraph::get.adjacency(yearly_graphLists[[1]][[l]],
                                                attr="weight" )}
  
  for(l in 1:L){
    nodeTensor_FC[[l]] <- igraph::get.adjacency(yearly_graphLists[[2]][[l]],
                                                attr="weight" )}
  
  
  #2. Built layer tensor
  layerTensor <- BuildLayersTensor(Layers=L, OmegaParameter=layerCouplingStrength, 
                                   MultisliceType=networkOfLayersType)
  
  
  
  #layerLabels <- c("Before", "short-term after", "long-term after")#1:L. Labels for each layer
  
  
  #3.1 Built the supra adjacency matrix (with different diagonals) of the multilayer network (unweighted)
  SupAdj_SC <- BuildSupraAdjacencyMatrixFromEdgeColoredMatrices(
    nodeTensor_SC, layerTensor, 
    L, Nodes)
  
  SupAdj_FC <- BuildSupraAdjacencyMatrixFromEdgeColoredMatrices(
    nodeTensor_FC, layerTensor, 
    L, Nodes)
  
  
  
  # Get multi centralities
  
  #strength................................................
  #
  list_multistrength_SC<- GetMultiStrength(SupAdj_SC,
                                           L,
                                           Nodes,
                                           isDirected=TRUE)%>%
    as.data.frame()
  
  
  list_multistrength_SC$names<-row.names(list_multistrength_SC)
  
  colnames(list_multistrength_SC)[1]<-"value"
  
  list_multistrength_SC$centrality<- "strength"  
  
  list_multistrength_SC$conn_type<- "SC"
  
  list_multistrength_SC$network_type<- "multilayer"
  
  
  
  
  #
  
  list_multistrength_FC<- GetMultiStrength(SupAdj_FC,
                                           L,
                                           Nodes,
                                           isDirected=TRUE)%>%
    as.data.frame()
  
  list_multistrength_FC$names<-row.names(list_multistrength_FC)
  
  colnames(list_multistrength_FC)[1]<-"value"
  
  list_multistrength_FC$centrality<- "strength"  
  
  list_multistrength_FC$conn_type<- "FC"
  
  list_multistrength_FC$network_type<- "multilayer"
  
  
  #PageRank............................................
  #                 
  
  list_multiPR_SC <- GetMultiPageRankCentrality(SupAdj_SC,L,Nodes)%>%
    as.data.frame()
  
  list_multiPR_SC$names<-row.names(list_multiPR_SC)
  
  colnames(list_multiPR_SC)[1]<-"value"
  
  list_multiPR_SC$centrality<- "PageRank"  
  
  list_multiPR_SC$conn_type<- "SC"
  
  list_multiPR_SC$network_type<- "multilayer"
  
  
  #Closeness............................................
  #                 
  
  list_multicloseness_SC <- GetMultiClosenessCentrality(SupAdj_SC,L,Nodes)%>%
    as.data.frame()
  
  list_multicloseness_SC$names<-row.names(list_multicloseness_SC)
  
  colnames(list_multicloseness_SC)[1]<-"value"
  
  list_multicloseness_SC$centrality<- "Closeness"  
  
  list_multicloseness_SC$conn_type<- "SC"
  
  list_multicloseness_SC$network_type<- "multilayer"
  
  #Hub............................................
  #                 
  
  list_multiHub_FC <- GetMultiHubCentrality(SupAdj_SC,L,Nodes)%>%
    as.data.frame()
  
  list_multiHub_FC$names<-row.names(list_multiHub_FC)
  
  colnames(list_multiHub_FC)[1]<-"value"
  
  list_multiHub_FC$centrality<- "Hub"  
  
  list_multiHub_FC$conn_type<- "FC"
  
  list_multiHub_FC$network_type<- "multilayer"
  
  
  multicentralities_SC_FC<- list_multistrength_SC %>% rbind(list_multistrength_FC,
                                                            list_multiPR_SC,list_multicloseness_SC,
                                                            list_multiHub_FC)
  
  
  
  return(multicentralities_SC_FC)
  
}
