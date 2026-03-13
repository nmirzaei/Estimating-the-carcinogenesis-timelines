# Estimating the carcinogenesis timelines
This repository contains MATLAB codes required to produce the results in "Estimating the carcinogenesis timelines in early-onset versus late-onset cancers and changes across birth cohorts"

# Using the codes
To produce the plots for each cancer type and each framework, use the Matlab files Plotting_{cancer type}_newdata.m

To produce the timelines for each cancer type and each framework, use the Matlab files Timeline_calc_{cancer type}_newdata.m

For identifiability and sensitivity analysis of the three mutation farmework, refer to the folder Three-mutation framework\Identifiability and Sensitivity. In this folder, you should go to the corresponding cancer folder the desired age cutoff, and then the desired cohort. For example, if you want to calculate the timelines for Breast cancer for the age cut-off <40 and the cohort of 1965-1969, you should go to:
Three-mutation framework\Identifiability and Sensitivity\Breast\less40\1965--1969
There, you run the file GA_cohorts.m

# Parameter estimation
The codes for parameter estimation are published with our earlier work. See: https://codeocean.com/capsule/6223220/tree/v1

This is a lengthy and computationally heavy process that needs to be done using a computing cluster 500 times. Therefore, we have provided the estimated parameters in this repository so that the codes can run on their own.

If you use our codes from the earlier study please cite the following paper:

Mohammad Mirzaei, Navid, et al. "Modeling early-onset cancer kinetics reveals changes in underlying risk and the impact of population screening." Cancer research 85.22 (2025): 4558-4570.
