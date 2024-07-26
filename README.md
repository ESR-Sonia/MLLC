# MLLC

Scripts for calculations of habitat connectivity based on the monolayer and multilayer centralities in a river-floodplain landscape.

## Functionality and background

This version of this script was written as part of the analysis performed in: Recinos S., Funk A., Tiwari S., Baldan D., Hein T. (2023). Multilayer Networks in Landscape Ecology: A Case Study to assess changes in aquatic habitat connectivity for flying and non-flying benthic Macroinvertebrates in a Danube floodplain [submitted for publication]. Institute of Hydrobiology and Aquatic Ecosystem Management, University of Natural Resources and Life Sciences, Vienna, Austria.

You may use this code to familiarize yourself with how to draw multilayer networks and calculate centralities using a list of edges and a list of nodes as input data. To follow the script “network_analysis.R” you will need to use the functions provided in the R folder.

## List of Files:

* R (folder)
  * get_monolayer_centralities.R
  * get_multicentralities.R
  * getlist_of_sequence_SC_FC_graphs.R
* Scripts (folder)
  * connected_components_under_daily_thresholds.R
  * mean_Closeness_under_daily_thresholds.R
  * mean_Hub_under_daily_thresholds.R
  * mean_PageRank_under_daily_thresholds.R
  * mean_weighted_degree_under_daily_thresholds.R
  * network_analysis.R
  * Input_data (folder)
    * edgeList_before.csv
    * edgeList_longT_after.csv
    * edgeList_shortT_after.csv
    * nodes_before.csv


## Acknowledgments

This work was funded by the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement number 859937 as part of the ‘i-CONN’ project. We also would like to thank Venetia Voutsa, Michalis Papadopoulos, Vicky Papadopoulou Lesta and Kay Schmidt for the conceptualization and development of earlier versions of the codes "connected_components_under_daily_thresholds.R" and the ones for mean centralities under daily thresholds.
