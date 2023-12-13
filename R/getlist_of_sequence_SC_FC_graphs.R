#' Creates a list of SC and FC graphs for the multilayer network
#'
#' @param edge_list object of class data.frame with the following columns: 
#' from, to, distance. This last column can be any cost/resistance edge attribute 
#' @param node_list object of class data.frame, with relevant node attributes as 
#' columns: Id, LHC.......
#' @param param_u numeric vector representing a upstream dispersal parameter. 
#' @param param_d numeric vector representing a downstream dispersal parameter. 
#' @details
#' \code{time_period = "before"} (default) is to be used to set maximum LHC 
#' threshold for the conditions before restoration 
#' \code{time_period = "after"} is to be used to set maximum LHC threshold for 
#' the conditions after restoration 
#' @references
#' Baldan, Damiano, et al. "Introducing ‘riverconn’: an R package to assess 
#' river connectivity indices." Environmental Modelling & Software 156 (2022).
#'   \href{https://doi.org/10.1016/j.envsoft.2022.105470}{doi:10.1016/j.envsoft.2022.105470}
#' @export
#' @references
#' Funk, A., Baldan, D., Bondar-Kunze, E., Brizuela, S. R., Kowal, J., & Hein, T. (2023). 
#' Connectivity as a driver of river-floodplain functioning: A dynamic, graph 
#' theoretic approach. Ecological Indicators, 154, 110877.
#' @export

getlist_of_sequence_SC_FC_graphs <- function(edge_list, node_list, param_u,param_d, 
                                             time_period ="before"){
  
  # get the vertices information
  graph <- graph_from_data_frame(select(edge_list, from, to), 
                                 #directed = TRUE,
                                 vertices = node_list)
  
  # add the length resistance/distance information to the edges information
  E(graph)$distance<-edge_list$distance
  
  
  graph_v_df <- igraph::as_data_frame(graph, what = "vertices") %>%
    dplyr::select(name, LHC)
  
  graph_e_df <- igraph::as_data_frame(graph, what = "edges") %>%
    dplyr::mutate(d_att = distance,
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
      mutate(flag_dir = "u"))
  
  
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
    dplyr::select(-c(disp_prob)) #this SC edge list has the same links as the FC, but it's weighted based on a measure of resistance
  
  
  
  ##Create yearly lists of graphs
  
  #1. Get the minimum LHC value between two pair of nodes
  #assuming that this link has at least this yearly LHC
  
  #unique(SC_edge_list$to)[!unique(SC_edge_list$to) %in% unique(graph_df_bidir_v$name)]#check if all nodes in edge file are in the node file
  
  SC_edge_list_2<-SC_edge_list %>%
    left_join(graph_df_bidir_v, by = c("from" = "name")) %>%
    mutate(hyd_up_val1_c=LHC) %>%
    left_join(graph_df_bidir_v, by = c("to" = "name")) %>%
    mutate(hyd_up_val2_c=LHC.y) %>%
    mutate(hc=pmin(hyd_up_val1_c,hyd_up_val2_c))%>% #min values of hc btw. each pair of nodes
    select(from,to,weight,value,cost,hc)
  
  FC_edge_list_2<-FC_edge_list %>%
    left_join(graph_df_bidir_v, by = c("from" = "name")) %>%
    mutate(hyd_up_val1_c=LHC) %>%
    left_join(graph_df_bidir_v, by = c("to" = "name")) %>%
    mutate(hyd_up_val2_c=LHC.y) %>%
    mutate(hc=pmin(hyd_up_val1_c,hyd_up_val2_c))%>% #min values of hc btw. each pair of nodes
    select(from,to,disp_prob,weight,hc)
  
  
  
  
  
  
  if(time_period == "before"){
    
    
    yearly_graphList_SC <- list()
    
    
    
    for(l in seq(1,171, by=10)){
      yearly_graphList_SC[[l]]<- graph_from_data_frame(SC_edge_list_2%>%
                                                         filter(hc >= l)%>%
                                                         mutate(weight=cost)%>%#i have to change the name of "cost" as "weight" so it can be recognized as a weighted graph in the nested for loop 
                                                         select(from, to,weight),
                                                       directed = TRUE,
                                                       vertices = graph_df_bidir_v)
      E(yearly_graphList_SC[[l]])$weight
      
      #return(yearly_graphList_SC)
      
    }
    
    yearly_graphList_FC <- list()
    
    
    
    
    for(l in seq(1,171, by=10)){
      
      yearly_graphList_FC[[l]]<- graph_from_data_frame(FC_edge_list_2%>%
                                                         filter(hc >= l)%>%
                                                         select(from, to,weight),
                                                       directed = F,
                                                       vertices = graph_df_bidir_v)
      E(yearly_graphList_FC[[l]])$weight
      
      #return(yearly_graphList_FC)
    }
  }
  
  
  
  
  
  if(time_period == "after"){
    
    
    yearly_graphList_SC <- list()
    
    
    for(l in seq(1,181, by=10)){
      yearly_graphList_SC[[l]]<- graph_from_data_frame(SC_edge_list_2%>%
                                                         filter(hc >= l)%>%
                                                         mutate(weight=cost)%>%#i have to change the name of "cost" as "weight" so it can be recognized as a weighted graph in the nested for loop 
                                                         select(from, to,weight),
                                                       directed = TRUE,
                                                       vertices = graph_df_bidir_v)
      E(yearly_graphList_SC[[l]])$weight
      
      #return(yearly_graphList_SC)
    }
    
    yearly_graphList_FC <- list()
    
    
    for(l in seq(1,181, by=10)){
      yearly_graphList_FC[[l]]<- graph_from_data_frame(FC_edge_list_2%>%
                                                         filter(hc >= l)%>%
                                                         select(from, to,weight),
                                                       directed = F,
                                                       vertices = graph_df_bidir_v)
      E(yearly_graphList_FC[[l]])$weight
      
      #return(yearly_graphList_FC)
    }
  }
  
  
  
  
  yearly_graphList_SC<- yearly_graphList_SC[lapply(yearly_graphList_SC,length)>0] 
  yearly_graphList_FC<- yearly_graphList_FC[lapply(yearly_graphList_FC,length)>0] 
  
  
  daily.graphs <- list(yearly_graphList_SC,yearly_graphList_FC)              # Store output in list
  return(daily.graphs)
  
  
}

