
#########################
##Script for network analysis in:
##'
##'Recinos S., Funk A., Tiwari S., Baldan D., Hein T. (2023). 
##'Multilayer Networks in Landscape Ecology: A
##'Case Study to assess changes in aquatic habitat
##'connectivity for flying and non-flying benthic
##'macroinvertebrates in a Danube floodplain
##'[submitted for publication]. Institute of Hydrobiology and 
##'Aquatic Ecosystem Management,
##'University of Natural Resources and Life Sciences, 
##'Vienna, Austria
##'
##'
##'Important references for network building:
##' @references
#' Baldan, Damiano, et al. "Introducing ‘riverconn’: an R package to assess 
#' river connectivity indices." Environmental Modelling & Software 156 (2022).
#'   \href{https://doi.org/10.1016/j.envsoft.2022.105470}{doi:10.1016/j.envsoft.2022.105470}
#' @export
#' 
##'Important references for calculating closeness:
##'@references
##'De Domenico, M., Porter, M. A., & Arenas, A. (2015). MuxViz: a tool for 
##'multilayer analysis and visualization of networks. Journal of Complex 
##'Networks, 3(2), 159-176.
##'  \href{https://doi.org/10.1093/comnet/cnu038}{}
##'@export


#########################



#--------------------------------------------------------------------------------------------------#
###
# Load required packages
list.of.packages <- c("igraph","tidyverse","data.table","muxViz","MASS","ggpubr")
invisible(lapply(list.of.packages, library, character.only = TRUE))


# Specify output directory, output_dir 

setwd()
#--------------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------------#
###Define sequential average closeness

getCloseness_of_dailyGraphs <- function(edge_list, node_list, 
                                        param_u,param_d,Nodes){
  
  #' @param edge_list object of class data.frame with the following columns: 
  #' from, to, distance. This last column can be any cost/resistance edge attribute 
  #' @param node_list object of class data.frame, with relevant node attributes as 
  #' columns: Id, LHC.......
  #' @param param_u numeric vector representing a upstream dispersal parameter. 
  #' @param param_d numeric vector representing a downstream dispersal parameter.
  #' @param Nodes numeric vector representing a downstream dispersal parameter.
  #' 
  #' get the vertices information
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
  
  ###FC networks
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
    dplyr::mutate(from=Var1, to=Var2,weight = (value/max(value)))%>%
    dplyr::select(from, to,value, weight)
  
  
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
  
  SC_edge_list_ <- FC_edge_list %>% select(from,to)%>%
    dplyr::left_join(SC_edge_list, by = c("from" = "Var1", "to" = "Var2")) %>%
    dplyr::filter(value!=Inf)%>%
    replace(is.na(.), 0)%>%
    dplyr::mutate(weight= (value/max(value))) %>%
    dplyr::select(from, to, value, weight) #this SC edge list has the same links as the FC, but it's weighted based on a measure of resistance
  
  ##
  ##Create yearly lists of graphs under daily thresholds
  ##
  
  #Get the minimum LHC value between two pair of nodes
  #assuming that an specific link has at least this LHC frequency
  
  
  SC_edge_list_2<-SC_edge_list_ %>%
    left_join(graph_df_bidir_v, by = c("from" = "name")) %>%
    mutate(hyd_up_val1_c=LHC) %>%
    left_join(graph_df_bidir_v, by = c("to" = "name")) %>%
    mutate(hyd_up_val2_c=LHC.y) %>%
    mutate(hc=pmin(hyd_up_val1_c,hyd_up_val2_c))%>% #min values of hc btw. each pair of nodes
    select(from,to,weight,value,hc)
  
  
  ##
  ##Calculate Closeness for a list of 365 networks
  ##
  
  ##with cost/resistance as link attribute
  
  
  list_closeness_SC <- list()
  SC_edge_list_3<-list()
  SC_graph_list<-list()
  DM<-list()
  
  for(i in 1:365){
    
    SC_edge_list_3[[i]]<-SC_edge_list_2%>%
      filter(hc > i)%>%select(-c(hc,value))
    
    SC_graph_list[[i]]<- igraph::graph_from_data_frame(SC_edge_list_3[[i]],
                                                       directed = TRUE,
                                                       vertices = graph_df_bidir_v)
    
    DM[[i]] <- igraph::distances(SC_graph_list[[i]], mode = "all")
    
    
    list_closeness_SC[[i]] <-
      unlist(lapply(1:Nodes, function(n)
        mean(1 / DM[[i]][n, ][-n])))%>%
      as.data.frame()
    
    list_closeness_SC[[i]]$names<-row.names(list_closeness_SC[[i]])
    colnames(list_closeness_SC[[i]])[1]<-"Closeness"
    
    list_closeness_SC[[i]]<- list_closeness_SC[[i]]%>% 
      summarise(mean_Closeness= mean(Closeness))
    
    list_closeness_SC[[i]]$threshold<- i
    
    
    
    
  }
  
  return(list(list_closeness_SC))
  
}


#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
#####import node files#####

nodes_before <- read.csv()
nodes_shortT_after <- read.csv()
nodes_longT_after <- read.csv() 



##import edge/link lists####

edgeList_before <- read.csv()
edgeList_shortT_after <- read.csv()
edgeList_longT_after <- read.csv()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
###calculate average Closeness for a list of 365 networks
param_u = 0 #for lowland passive aquatic dispersers
param_d = 3
Nodes= 886


############
###Before
comp1<-getCloseness_of_dailyGraphs(edge_list=edgeList_before, 
                                   node_list=nodes_before,
                                   param_u,param_d,Nodes)

comp1_SC<-data.frame(matrix(unlist(comp1[[1]]), nrow=length(comp1[[1]]), byrow=TRUE)) %>%
  mutate(mean_Closeness= X1, threshold= X2, Period= "Before")%>% select(mean_Closeness, threshold, Period)%>%
  group_by(threshold)


###short-term after,
comp2<-getCloseness_of_dailyGraphs(edge_list=edgeList_shortT_after, 
                                   node_list=nodes_shortT_after,
                                   param_u,param_d,Nodes)

comp2_SC<-data.frame(matrix(unlist(comp2[[1]]), nrow=length(comp2[[1]]), byrow=TRUE)) %>%
  mutate(mean_Closeness= X1, threshold= X2, Period= "Short-term after")%>%select(mean_Closeness, threshold, Period)%>%
  group_by(threshold)


###Long-term after,
comp3<-getCloseness_of_dailyGraphs(edge_list=edgeList_longT_after, 
                                   node_list=nodes_longT_after,
                                   param_u,param_d,Nodes)

comp3_SC<-data.frame(matrix(unlist(comp3[[1]]), nrow=length(comp3[[1]]), byrow=TRUE)) %>%
  mutate(mean_Closeness= X1, threshold= X2, Period= "Long-term after")%>%select(mean_Closeness, threshold, Period)%>%
  group_by(threshold)




##thresholds

#for before
#seq(1,171, by=10)
#thresholds: 1  11  21  31  41  51  61  71  81  91 101 111 121 131 141 151 161 171



#for after
#seq(1,181, by=10)
#thresholds: 1  11  21  31  41  51  61  71  81  91 101 111 121 131 141 151 161 171 181
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
###Plot average Closeness for a list of 365 networks


###Plots
p<-ggplot()+
  geom_step(data=comp1_SC, mapping=aes(x = threshold, y= mean_Closeness, color= Period),
            direction = "vh",lwd = 1.2 )+
  geom_step(data=comp2_SC, mapping=aes(x = threshold, y= mean_Closeness,color=Period),
            direction = "vh",lwd = 1.2   )+
  geom_step(data=comp3_SC, mapping=aes(x = threshold, y= mean_Closeness,color=Period),
            direction = "vh",lwd = 1.2  )+
  xlab("Minimum days of LHC") +
  ylab("Mean Closeness") +
  #ggtitle("Normalized dPC index for different floodplain regions")+
  #xlim(0, 363)+
  geom_vline(xintercept=181,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=171,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=161,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=151,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=141,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=131,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=121,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=111,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=101,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=91,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=81,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=71,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=61,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=51,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=41,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=31,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=21,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=11,linetype='dashed', color="#972D15", size=0.5)+
  geom_vline(xintercept=1.5,linetype='dashed', color="#972D15", size=0.5)+
  guides(colour=guide_legend(ncol = 2, title.position="top"))+
  scale_x_continuous(trans = scales::reverse_trans(), 
                     expand = c(0, 0),limits = c(200,1))+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 0.15))+
  expand_limits(x = 0, y = 0)+
  theme_bw()+
  theme(legend.position="top",axis.text.x = element_text(size = 12, angle = 45, hjust=1)) +
  scale_color_manual(values= c("#EBCC2A", "#78B7C5", "#3B9AB2"))+
  ggtitle(expression("Mean closeness at daily LHC thresholds")) 

p

png(file = "mean_closeness_365thresholds.png", res = 300, width = 2000, height = 1315)
plot (p)
dev.off()



