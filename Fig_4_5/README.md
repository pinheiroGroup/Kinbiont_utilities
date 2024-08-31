
In this folder it is possible to find the scripts and data to replicate the analysis of decision tree of data from "High-throughput characterization of bacterial responses to complex mixtures of chemical pollutants", 2024, Smith et al, DOI: 10.1038/s41564-024-01626-9


Please note that to perform the fit it is necessary to extract the data with cleaned data.
To match isolate number and names, please look at the original github https://github.com/smithtp/isolate-chem-mixtures/


1. 'loop_chem_isolates_analysis_NL.jl' performs the NL fit with Richards model of all experiments
2. 'Loop_DT_with_maxdeph.jl ' performs the decision tree regression on all strains at different depth. It returns the Impurity scores, the Impurity rank, and the 10-fold cross validation R^2 for each parameter depth and strain
3. 'Loop_DT_with_maxdeph.jl ' performs the decision tree regression on mixture and N. soli at different depth. It returns the impurity scores, the Impurity rank, and the 10-flod cross validation R^2 for each parameter depth and the plots of the trees.

