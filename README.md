# MLLC

Scripts for calculations of habitat connectivity based on the monolayer and multilayer centralities in a river-floodplain landscape.

## Functinality and background

The first version of this script was written as part of the analysis performed in: Recinos S., Funk A., Tiwari S., Baldan D., Hein T. (2023). Multilayer Networks in Landscape Ecology: A Case Study to assess changes in aquatic habitat connectivity for flying and non-flying benthic macroinvertebrates in a Danube floodplain [submitted for publication]. Institute of Hydrobiology and Aquatic Ecosystem Management, University of Natural Resources and Life Sciences, Vienna, Austria.

You may use this code to familiarize yourself with how to draw multilayer networks and calculate centralities using a list of edges and a list of nodes as input data. To follow the script “network_analysis.R” you will need to use the functions provided in the R folder.

## List of Files:

1.	Scripts (folder)
1.1.	connected_components_under_daily_thresholds.R
1.2.	mean_Closeness_under_daily_thresholds.R
1.3.	mean_Hub_under_daily_thresholds.R
1.4.	mean_PageRank_under_daily_thresholds.R
1.5.	mean_weighted_degree_under_daily_thresholds.R
1.6.	network_analysis.R
1.7.	Input_data (folder)
1.7.1.	edgeList_before.csv
1.7.2.	edgeList_longT_after.csv
1.7.3.	edgeList_shortT_after.csv
1.7.4.	nodes_before.csv
1.7.5.	nodes_longT_after.csv
1.7.6.	nodes_shortT_after.csv
2.	R (folder)
2.1.	get_monolayer_centralities.R
2.2.	get_multicentralities.R
2.3.	getlist_of_sequence_SC_FC_graphs.R

## Acknowledgments

This work was funded by the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement number 859937 as part of the ‘i-CONN’ project. We also would like to thank Venetia Voutsa, Michalis Papadopoulos, Vicky Papadopoulou Lesta and Kay Schmidt for the conceptualization and development of earlier versions of the codes "connected_components_under_daily_thresholds.R" and the ones for mean centralities under daily thresholds.
