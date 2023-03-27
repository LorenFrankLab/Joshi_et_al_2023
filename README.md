# Joshi_et_al_2023
Analysis code for Joshi et al 2023 Dynamic synchronization between hippocampal representations and stepping.

	This repo contains code to analyze and generate figures for Joshi et al., 2023 "Dynamic synchronization between hippocampal representations and stepping.". 
    Data was acquired using Trodes1.8.2 and 30-tet DKR drive versions for five rats, extracted using Trodes export functions, then processed and analyzed in MATLAB 2020a.
    Contact Abhilasha Joshi (abhilasha.joshi@ucsf.edu) with any questions. 

Repos used to preprocess and analyze the data:

	Data processing:  https://bitbucket.org/franklab/trodes2ff_shared/ 
	Analysis: https://bitbucket.org/franklab/filterframework_shared/, https://github.com/tjd2002/tjd-shared-code and https://github.com/edeno/pose_analysis. 
	Figures: were saved as pdf or eps and formatted in Adobe Illustrator. 
	Data in NWB format is available on the DANDI: https://dandiarchive.org/dandiset/000410/draft/. Mat files to generate figures are on https://zenodo.org/deposit/7615939. 
 
Figure 1

	- plot_example_lineartrack.m
	- dfa_clusterfilter
	- Figure1ACD.mat

	- dfs_plotpower_outbound_inbound_peranimal.m
	- dfa_steptheta_wtrack.m
	- smooth_dlc_stepcycles_wtrack.m
	- Figure1ACD.mat

	- dfs_speed_steps_theta.m
	- dfa_speed_steps_theta_forelimb
	- Figure1EF_EDFigure2.mat

	- dfs_speed_steps_theta_perepoch
	- Figure1EF_EDFigure2.mat
	
Figures 2 and 3

	- dfs_reldist_analysis_allcombination.m
	- dfa_ahead_behind_distance_muareldist.m
	- dfa_ahead_behind_distance_muareldist_lintrack.m
	- Figure23_plant_triggered_avgs.mat
	- Figure2D_EDFigure7BC_decode_to_animal_dist.mat

ED Figure 1

	- validate_plant_times.m
	- EDFigure1.mat


ED Figure 2

	- dfs_speed_steps_theta.m
	- dfa_speed_steps_theta_forelimb
	- Figure1EF_EDFigure2.mat

	- dfs_speed_steps_theta_perepoch
	- Figure1EF_EDFigure2.mat

ED Figure 3

	- dfs_gait_motifs_analysis.m
	- EDFigure3.mat
	- dfa_steptheta_wtrack_gait.m

ED Figures 4 and 5

	- dfs_reldist_analysis_allcombination.m
	- dfs_mua_analysis_allcombination.m
	- EDFigure45C_EDFigure7DE_mua_mod_score.mat
	- dfa_ahead_behind_distance_muareldist.m
	- dfa_ahead_behind_distance_muareldist_lintrack.m

ED Figure 6

	- dfs_thetaseq_length.m
	- EDFigure6.mat

ED Figure 7

	- dfs_reldist_analysis_allcombination.m
	- dfs_mua_analysis_allcombination.m
	- Figure2D_EDFigure7BC_decode_to_animal_dist.mat
	- EDFigure45C_EDFigure7DE_mua_mod_score.mat
