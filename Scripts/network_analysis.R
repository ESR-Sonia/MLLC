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
##'Important references for calculating multilayer centralities:
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
list.of.packages <- c("igraph","tidyverse","data.table","muxViz","MASS","heatmaply")
invisible(lapply(list.of.packages, library, character.only = TRUE))


# Specify output directory, output_dir 

setwd()
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
###import node files#####

nodes_before <- read.csv()
nodes_shortT_after <- read.csv()
nodes_longT_after <- read.csv()



##import edge/link lists####

edgeList_before <- read.csv()
edgeList_shortT_after <- read.csv()
edgeList_longT_after <- read.csv()

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
#####Get monolayer centralities#########


param_u = 0 #upstream dispersal parameter for lowland passive aquatic dispersers
param_d = 3 #downstream dispersal parameter for lowland passive aquatic dispersers
Nodes= 886 #number of nodes

Distance_monocentralities_96<-get_monolayer_centralities(edge_list=edgeList_before, 
                                                         node_list=nodes_before, 
                                                         param_u,param_d,Nodes)

Distance_monocentralities_96<-Distance_monocentralities_96%>%
  select(names, everything())%>% mutate(year="1996")


#separate for table 
#
Distance_strength1<- Distance_monocentralities_96 %>%
  filter(centrality=="strength", conn_type=="SC", network_type=="monolayer")%>%
  mutate(WD_SC96_netdist=value)%>%select(names,WD_SC96_netdist)

Distance_strength<- Distance_monocentralities_96 %>%
  filter(centrality=="strength", conn_type=="FC", network_type=="monolayer")%>%
  mutate(WD_FC96_netdist=value)%>%select(names,WD_FC96_netdist)

Distance_PR1<- Distance_monocentralities_96 %>%
  filter(centrality=="PageRank")%>%
  mutate(PR96_netdist=value)%>%select(names,PR96_netdist)

Distance_closeness1<- Distance_monocentralities_96 %>%
  filter(centrality=="Closeness")%>%
  mutate(closeness96_netdist=value)%>%select(names,closeness96_netdist)

Distance_Hub<- Distance_monocentralities_96 %>%
  filter(centrality=="Hub")%>%
  mutate(Hub96_netdist=value)%>%select(names,Hub96_netdist)

Distance_monocentralities96_PLS_plots<-list(Distance_strength1,Distance_strength,Distance_PR1,
                                            Distance_closeness1,Distance_Hub)

Distance_monocentralities96_PLS_plots<-Distance_monocentralities96_PLS_plots %>% 
  reduce(full_join, by='names')



#1999

Distance_monocentralities_99<-get_monolayer_centralities(edgeList_shortT_after,
                                                         nodes_shortT_after, 
                                                         param_u,param_d,
                                                         Nodes)

Distance_monocentralities_99<-Distance_monocentralities_99%>%
  select(names, everything())%>% mutate(year="1999")


#separate for table 
#
Distance_strength1<- Distance_monocentralities_99 %>%
  filter(centrality=="strength", conn_type=="SC", network_type=="monolayer")%>%
  mutate(WD_SC99_netdist=value)%>%select(names,WD_SC99_netdist)

Distance_strength2<- Distance_monocentralities_99 %>%
  filter(centrality=="strength", conn_type=="FC", network_type=="monolayer")%>%
  mutate(WD_FC99_netdist=value)%>%select(names,WD_FC99_netdist)

Distance_PR2<- Distance_monocentralities_99 %>%
  filter(centrality=="PageRank")%>%
  mutate(PR99_netdist=value)%>%select(names,PR99_netdist)

Distance_closeness2<- Distance_monocentralities_99 %>%
  filter(centrality=="Closeness")%>%
  mutate(closeness99_netdist=value)%>%select(names,closeness99_netdist)

Distance_Hub<- Distance_monocentralities_99 %>%
  filter(centrality=="Hub")%>%
  mutate(Hub99_netdist=value)%>%select(names,Hub99_netdist)


Distance_monocentralities99_PLS_plots<-list(Distance_strength1,Distance_strength2,Distance_PR2,
                                            Distance_closeness2,Distance_Hub)

Distance_monocentralities99_PLS_plots<-Distance_monocentralities99_PLS_plots %>% 
  reduce(full_join, by='names')




#2020


Distance_monocentralities_2020<-get_monolayer_centralities(edgeList_longT_after, 
                                                           nodes_longT_after, 
                                                           param_u,param_d,
                                                           Nodes)

Distance_monocentralities_2020<-Distance_monocentralities_2020%>%
  select(names, everything())%>% mutate(year="2020")


#separate for table 
#
Distance_strength1<- Distance_monocentralities_2020 %>%
  filter(centrality=="strength", conn_type=="SC", network_type=="monolayer")%>%
  mutate(WD_SC2020_netdist=value)%>%select(names,WD_SC2020_netdist)

Distance_strength3<- Distance_monocentralities_2020 %>%
  filter(centrality=="strength", conn_type=="FC", network_type=="monolayer")%>%
  mutate(WD_FC2020_netdist=value)%>%select(names,WD_FC2020_netdist)

Distance_PR3<- Distance_monocentralities_2020 %>%
  filter(centrality=="PageRank")%>%
  mutate(PR2020_netdist=value)%>%select(names,PR2020_netdist)

Distance_closeness3<- Distance_monocentralities_2020 %>%
  filter(centrality=="Closeness")%>%
  mutate(closeness2020_netdist=value)%>%select(names,closeness2020_netdist)

Distance_Hub<- Distance_monocentralities_2020 %>%
  filter(centrality=="Hub")%>%
  mutate(Hub2020_netdist=value)%>%select(names,Hub2020_netdist)


Distance_monocentralities2020_PLS_plots<-list(Distance_strength1,Distance_strength3,Distance_PR3,
                                              Distance_closeness3,Distance_Hub)

Distance_monocentralities2020_PLS_plots<-Distance_monocentralities2020_PLS_plots %>% 
  reduce(full_join, by='names')


#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#

##### Multilayer centralities ###


#####Step 1.: Creating list of graph objects daily LHC thresholds (layers)#####

#directed distances
Distance_graphs_96<-getlist_of_sequence_SC_FC_graphs(edgeList_before, nodes_before, param_u,param_d, 
                                                     time_period ="before")
Distance_graphs_99<-getlist_of_sequence_SC_FC_graphs(edgeList_shortT_after, nodes_shortT_after, 
                                                     param_u,param_d, 
                                                       time_period ="after")
Distance_graphs_2020<-getlist_of_sequence_SC_FC_graphs(edgeList_longT_after, nodes_longT_after, 
                                                       param_u,param_d, 
                                                       time_period ="after")



#####Step .2: Building Multilayer networks and calculating multilayer centralities#####



L_1<- 18 #number of layers for "before"
L_2<- 19 #number of layers for "after"
#Nodes<- 886 #number of nodes
layerCouplingStrength <- 1 #This is an unweighted multilayer (with weighted intralayer links but unweighted interlayer links)
networkOfLayersType <- "categorical"



##network distances............................

#96
Distance_multicentralities_96<-get_multicentralities(Distance_graphs_96,
                                                     L_1,Nodes,
                                                     layerCouplingStrength,
                                                     networkOfLayersType)

Distance_multicentralities_96<-Distance_multicentralities_96%>%
  select(names, everything())%>% mutate(year="1996")


#separate for table 
#
Distance_strength1<- Distance_multicentralities_96 %>%
  filter(centrality=="strength", conn_type=="SC")%>%
  mutate(multiWD_SC96_netdist=value)%>%select(names,multiWD_SC96_netdist)

Distance_strength<- Distance_multicentralities_96 %>%
  filter(centrality=="strength", conn_type=="FC")%>%
  mutate(multiWD_FC96_netdist=value)%>%select(names,multiWD_FC96_netdist)

Distance_PR<- Distance_multicentralities_96 %>%
  filter(centrality=="PageRank")%>%
  mutate(multiPR96_netdist=value)%>%select(names,multiPR96_netdist)

Distance_closeness<- Distance_multicentralities_96 %>%
  filter(centrality=="Closeness")%>%
  mutate(multicloseness96_netdist=value)%>%select(names,multicloseness96_netdist)

Distance_Hub<- Distance_multicentralities_96 %>%
  filter(centrality=="Hub")%>%
  mutate(multiHub96_netdist=value)%>%select(names,multiHub96_netdist)

Distance_multicentralities96_PLS_plots<-list(Distance_strength1,Distance_strength,Distance_PR,
                                             Distance_closeness,Distance_Hub)

Distance_multicentralities96_PLS_plots<-Distance_multicentralities96_PLS_plots %>% 
  reduce(full_join, by='names')


#99
Distance_multicentralities_99<-get_multicentralities(Distance_graphs_99,
                                                     L=L_2,Nodes,
                                                     layerCouplingStrength,
                                                     networkOfLayersType)

Distance_multicentralities_99<-Distance_multicentralities_99%>%
  select(names, everything())%>% mutate(year="1999")


#separate for table 
#
Distance_strength1<- Distance_multicentralities_99 %>%
  filter(centrality=="strength", conn_type=="SC")%>%
  mutate(multiWD_SC99_netdist=value)%>%select(names,multiWD_SC99_netdist)

Distance_strength<- Distance_multicentralities_99 %>%
  filter(centrality=="strength", conn_type=="FC")%>%
  mutate(multiWD_FC99_netdist=value)%>%select(names,multiWD_FC99_netdist)

Distance_PR<- Distance_multicentralities_99 %>%
  filter(centrality=="PageRank")%>%
  mutate(multiPR99_netdist=value)%>%select(names,multiPR99_netdist)

Distance_closeness<- Distance_multicentralities_99 %>%
  filter(centrality=="Closeness")%>%
  mutate(multicloseness99_netdist=value)%>%select(names,multicloseness99_netdist)

Distance_Hub<- Distance_multicentralities_99 %>%
  filter(centrality=="Hub")%>%
  mutate(multiHub99_netdist=value)%>%select(names,multiHub99_netdist)

Distance_multicentralities99_PLS_plots<-list(Distance_strength1,Distance_strength,Distance_PR,
                                             Distance_closeness,Distance_Hub)

Distance_multicentralities99_PLS_plots<-Distance_multicentralities99_PLS_plots %>% 
  reduce(full_join, by='names')




#2020
Distance_multicentralities_2020<-get_multicentralities(Distance_graphs_2020,
                                                       L=L_2,Nodes,
                                                       layerCouplingStrength,
                                                       networkOfLayersType)

Distance_multicentralities_2020<-Distance_multicentralities_2020%>%
  select(names, everything())%>% mutate(year="2020")


#separate for table 
#
Distance_strength1<- Distance_multicentralities_2020 %>%
  filter(centrality=="strength", conn_type=="SC")%>%
  mutate(multiWD_SC2020_netdist=value)%>%select(names,multiWD_SC2020_netdist)

Distance_strength<- Distance_multicentralities_2020 %>%
  filter(centrality=="strength", conn_type=="FC")%>%
  mutate(multiWD_FC2020_netdist=value)%>%select(names,multiWD_FC2020_netdist)

Distance_PR<- Distance_multicentralities_2020 %>%
  filter(centrality=="PageRank")%>%
  mutate(multiPR2020_netdist=value)%>%select(names,multiPR2020_netdist)

Distance_closeness<- Distance_multicentralities_2020 %>%
  filter(centrality=="Closeness")%>%
  mutate(multicloseness2020_netdist=value)%>%select(names,multicloseness2020_netdist)

Distance_Hub<- Distance_multicentralities_2020 %>%
  filter(centrality=="Hub")%>%
  mutate(multiHub2020_netdist=value)%>%select(names,multiHub2020_netdist)

Distance_multicentralities2020_PLS_plots<-list(Distance_strength1,Distance_strength,Distance_PR, 
                                               Distance_closeness,Distance_Hub)

Distance_multicentralities2020_PLS_plots<-Distance_multicentralities2020_PLS_plots %>% 
  reduce(full_join, by='names')





####create long table####


mono_multi_centralities_long<-rbind(Distance_monocentralities_96,
                                        Distance_monocentralities_99,
                                        Distance_monocentralities_2020,
                                        Distance_multicentralities_96,
                                        Distance_multicentralities_99,
                                        Distance_multicentralities_2020
                                        )

mono_multi_centralities_long<-mono_multi_centralities_long%>% 
  select(names,centrality,conn_type,network_type,value)%>%
  write_excel_csv("mono_multi_centralities_long.csv", delim = ",")

