##'Network construction details:
##'
##'In monolayer networks links consists in possible dispersal routes (OD-links) 
##'for a passive lowland dispersers, following the approach of:
##'
##'Baldan, D., Cunillera-Montcusí, D., Funk, A., & Hein, T. (2022). Introducing 
##'‘riverconn’: an R package to assess river connectivity indices. Environmental 
##'Modelling & Software, 156, 105470.
##'
##'In this approach, links between direct neighbouring and non-direct 
##'neighbouring nodes are estimated based on the likelihood 
##'for dispersal which decays exponentially with distance/resistance 
##'(see: https://doi.org/10.1016/j.ecolmodel.2010.06.017).
##'
##'##'The resulting networks were weighted based on:
##'The probability for dispersal are set as FC weights, which are given by the 
##'exponential kernel.
##'Network (in-channel) distances are used as cost attribute to weight the final 
##'OD-links.
##' @param edge_list dataframe specifying the ID of the starting and end node 
##' (from, to), with cost attribute as the column "distance" 
##' @param node_list dataframe specifying the ID of the nodes and yearly 
##' frequencies of inundation (LHC column)
##' @param param_u numeric vector with upstream dispersal parameter for passive 
##' aquatic dispersers
##' @param param_d numeric vector with downstream dispersal parameter for 
##' passive aquatic dispersers
##' @param Nodes numeric vector specifying the number of nodes in the node list 

get_monolayer_centralities<-function(edge_list, node_list, param_u,param_d,Nodes){
  
  # get the vertices information
  graph <- graph_from_data_frame(dplyr::select(edge_list, from, to), 
                                 #directed = TRUE,
                                 vertices = node_list)
  # add the length resistance/distance information to the edges information
  E(graph)$distance<-edge_list$distance
  
  
  graph_v_df <- igraph::as_data_frame(graph, what = "vertices") %>%
    dplyr::select(name, LHC)
  
  
  
  graph_e_df <- igraph::as_data_frame(graph, what = "edges") %>%
    dplyr::mutate(d_att = edge_list$distance,
                  flag_dir = "n")
  
  
  # create intermediate graph that has the distances
  graph_intermediate <- igraph::graph_from_data_frame(d = graph_e_df, vertices = graph_v_df)
  
  # If directionality is "asymmetric" the graph is made undirectional
  # Bidirectional edge list
  graph_df_bidir<- rbind(
    # "downstream" graph: graph with directions going downstream, uses pass_d as attribute
    igraph::as_data_frame(graph_intermediate, what = "edges") %>%
      dplyr::mutate(flag_dir = "d"),
    # "upstream" graph: graph with directions going upstream, uses pass_u as attribute
    igraph::as_data_frame(graph_intermediate, what = "edges") %>%
      dplyr::rename(from = to, to = from) %>%
      mutate(flag_dir = "u")
  )
  
  
  # Vertices dataframe
  
  graph_df_bidir_v <- igraph::as_data_frame(graph, what = "vertices") %>%
    dplyr::relocate(name)
  
  
  # Create the bidirectional graph
  
  graph_output <- igraph::graph_from_data_frame(
    d = graph_df_bidir,
    vertices = graph_df_bidir_v )
  
  
  
  #after creating bidirectional graph get the shortest paths separately
  
  
  #
  # asymmetric dispersal: I use both distances
  #
  # Create graph for upstream movement
  
  graph_dodgr_u <- igraph::as_data_frame(graph_output, what = "edges") %>%
    filter(flag_dir == "u") %>%
    rename(dist = d_att) %>%
    select(from, to, dist)
  
  
  #make directed edge list (add zeros to upper triangular)
  edgelist_u <- rbind(graph_dodgr_u,
                      graph_dodgr_u %>%
                        rename(from = to, to = from) %>%
                        mutate(dist = 0))
  
  
  
  graph_dodgr_u <- igraph::graph_from_data_frame(
    d = edgelist_u,
    vertices = graph_df_bidir_v )
  
  
  # Calculate all shortest paths for upstream movement
  
  
  Bij_mat_u <- igraph::distances(
    graph_dodgr_u,
    v = V(graph_dodgr_u),
    to = V(graph_dodgr_u),
    mode = c("out"),
    weights = E(graph_dodgr_u)$dist,
    algorithm = c("dijkstra")
  )
  
  
  # Create dodgr graph for downstream movement
  
  graph_dodgr_d <- igraph::as_data_frame(graph_output, what = "edges") %>%
    filter(flag_dir == "d") %>%
    rename(dist = d_att) %>%
    select(from, to, dist)
  
  #make directed edge list (add zeros to lower triangular)
  edgelist_d <- rbind(graph_dodgr_d,
                      graph_dodgr_d %>%
                        rename(from = to, to = from) %>%
                        mutate(dist = 0))
  
  graph_dodgr_d <- igraph::graph_from_data_frame(
    d = edgelist_d,
    vertices = graph_df_bidir_v )
  
  
  # Calculate all shortest paths for downstream movement
  
  Bij_mat_d <- igraph::distances(
    graph_dodgr_d,
    v = V(graph_dodgr_d),
    to = V(graph_dodgr_d),
    mode = c("out"),
    weights = E(graph_dodgr_d)$dist,
    algorithm = c("dijkstra")
  )
  
  # post process FC matrix if exponential decay
  
  Bij_mat = (param_u^Bij_mat_u) * (param_d^Bij_mat_d )
  
  
  FC_edge_list <- cbind(expand.grid(dimnames(Bij_mat)), value = as.vector(Bij_mat))
  
  
  FC_edge_list <- FC_edge_list %>%
    dplyr::filter(value>0,value!=Inf)%>%
    na.omit()%>%
    dplyr::mutate(from=Var1, to=Var2,disp_prob=value, weight = (value/max(value)))%>%
    dplyr::select(from, to,disp_prob, weight)
  
  
  ##Create FC network with weights: dispersal probabilities
  graph_FC<- graph_from_data_frame(select(FC_edge_list, from, to,disp_prob), 
                                   directed = F,
                                   vertices = graph_df_bidir_v)
  
  E(graph_FC)$weight<-FC_edge_list$weight
  
  
  ##create directed graph with same number of links, but different weights
  
  
  ##SC weights: resistances/distances
  
  ## calculate shortest paths with both upper and lower triangular  
  
  Bij_mat_all <- igraph::distances(
    graph_output,
    v = V(graph_output),
    to = V(graph_output),
    mode = c("all"),
    weights = E(graph_output)$d_att,
    algorithm = c("dijkstra")
  )
  
  
  SC_edge_list <- cbind(expand.grid(dimnames(Bij_mat_all)), value = as.vector(Bij_mat_all))
  
  
  SC_edge_list <- FC_edge_list %>%
    dplyr::left_join(SC_edge_list, by = c("from" = "Var1", "to" = "Var2")) %>%
    dplyr::mutate(cost= (value/max(value))) %>% 
    dplyr::select(-c(disp_prob))
  
  
  graph_SC<- graph_from_data_frame(select(SC_edge_list, from, to,cost), 
                                   directed = TRUE,
                                   vertices = graph_df_bidir_v)
  
  E(graph_SC)$cost<-SC_edge_list$cost
  
  
  
  ##Calculate centralities
  
  
  
  #strength............................................
  
  #SC
  
  #Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro 
  #Vespignani: The architecture of complex weighted networks, Proc. Natl. Acad. 
  #Sci. USA 101, 3747 (2004)
  
  #
  
  list_strength_SC <- igraph::strength(graph_SC,mode = c("all"))%>%
    as.data.frame()
  
  list_strength_SC$names<-row.names(list_strength_SC)
  
  colnames(list_strength_SC)[1]<-"value"
  
  list_strength_SC$centrality<- "strength"  
  
  list_strength_SC$conn_type<- "SC"
  
  list_strength_SC$network_type<- "monolayer"
  
  
  #FC
  
  
  #
  
  list_strength_FC <- igraph::strength(graph_FC,mode = c("all"))%>%
    as.data.frame()
  
  list_strength_FC$names<-row.names(list_strength_FC)
  
  colnames(list_strength_FC)[1]<-"value"
  
  list_strength_FC$centrality<- "strength"  
  
  list_strength_FC$conn_type<- "FC"
  
  list_strength_FC$network_type<- "monolayer"
  
  
  
  #PageRank............................................
  
  #Sergey Brin and Larry Page: The Anatomy of a Large-Scale Hypertextual Web 
  #Search Engine. Proceedings of the 7th World-Wide Web Conference, Brisbane, 
  #Australia, April 1998.
  
  #SC
  
  #
  
  list_PR_SC <- igraph::page_rank(graph_SC,vids=V(graph_SC),
                                  directed = T)$vector%>%
    as.data.frame()
  
  list_PR_SC$names<-row.names(list_PR_SC)
  
  colnames(list_PR_SC)[1]<-"value"
  
  list_PR_SC$centrality<- "PageRank"  
  
  list_PR_SC$conn_type<- "SC"
  
  list_PR_SC$network_type<- "monolayer"
  
  
  
  #closeness............................................
  
  #Opsahl, T., Agneessens, F., Skvoretz, J. (2010). Node centrality in weighted 
  #networks: Generalizing degree and shortest paths. Social Networks 32, 245-251
  #https://toreopsahl.com/2010/03/20/closeness-centrality-in-networks-with-disconnected-components/
  
  #SC
  
  #
  DM <- igraph::distances(graph_SC, mode = "all")
  
  list_closeness_SC <-
    unlist(lapply(1:Nodes, function(n)
      mean(1 / DM[n, ][-n])))%>%
    as.data.frame()
  
  #0 if disconnected, 1 if connected to all other nodes
  
  list_closeness_SC$names<-row.names(list_closeness_SC)
  
  colnames(list_closeness_SC)[1]<-"value"
  
  list_closeness_SC$centrality<- "Closeness"  
  
  list_closeness_SC$conn_type<- "SC"
  
  list_closeness_SC$network_type<- "monolayer"
  
  
  
  #Hub.........................................
  
  #J. Kleinberg. Authoritative sources in a hyperlinked environment. Proc. 9th 
  #ACM-SIAM Symposium on Discrete Algorithms, 1998. Extended version in Journal 
  #of the ACM 46(1999). Also appears as IBM Research Report RJ 10076, May 1997.
  
  
  #FC
  
  
  #
  
  list_hub_FC <- igraph::hub_score(graph_FC, scale = TRUE)$vector%>%
    as.data.frame()
  
  list_hub_FC$names<-row.names(list_hub_FC)
  
  colnames(list_hub_FC)[1]<-"value"
  
  list_hub_FC$centrality<- "Hub"  
  
  list_hub_FC$conn_type<- "FC"
  
  list_hub_FC$network_type<- "monolayer"
  
  
  #Merge in a single dataframe
  monocentralities_SC_FC<- list_strength_SC %>% rbind(list_strength_FC,
                                                      list_PR_SC,
                                                      list_closeness_SC,
                                                      list_hub_FC)
  
  
  
}
