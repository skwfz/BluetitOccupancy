# BluetitOccupancy

This repository contains R and Stan code for a project I did as a part of the Bayesian Data Analysis Course at Aalto University. The goal was to analyze the occupancy of the Eurasian Blue Tit in Finland using fully Bayesian methods.

The code is not necessarily very well organized, but occupancy_data_preparation.R and land_covariates.R contain code for data preprocessing, which were largely based on the tutorial at https://cornelllabofornithology.github.io/ebird-best-practices/. The file runnin_models.R contains code for running Stan models and plotting some results. bird_distribution.R plots an occupancy map using the Stan models.

Note that all of the data is missing here, and they are available from Ebird (https://ebird.org/) on request. I got the satellite data from NASA websites (following the tutorial). 
