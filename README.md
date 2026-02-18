# smurf_stat
This repository contains the code for the paper *Acute Smurf mortality and inter-phase dependence in Drosophila and mice identified through comprehensive modelling and statistical analysis of two-phase ageing* by L.Breuil, M.Doumic, S.Kaakaï and M.Rera

The files Estimation.R and Simulation.R contain the source R code and functions. 

The file Smurf_stat_study.ipynb contains the code generating all of the parameter values and generates the figures of the article. 

The Figures of the article are in the folder Figure. 

The data used in the article is in the folder Data : times_exp is the main Drosophila data on which the the study was conducted, C57males, C57femelles and AJRJ contain the mouse data, surv_lines contains the survival data for DGRP line 377 
and Var_smurf_mx, Var_death_tot and Var_death_mx contain computed variances for confidence intervals which are very long to compute. 