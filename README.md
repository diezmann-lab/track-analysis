A repository for code used to analyze single-molecule tracking data. This code goes from a .csv file including XYZT data, is refined using swift v. 0.4.3, and is then analyzed in MATLAB.

Pipeline:
1. Export data
Export single-molecule localization data.
Inputs: Experimental sample.
Outputs: .csv file of localizations including x,y positions and acquisition frame.
Dependencies: none

2. Assign trajectories
Generate trajectories from localizations using Swift. Calculate start times and duration.
Inputs: .csv file with localizations.
Outputs: .csv file with trajectories.
Dependencies: none

3. Prefilter trajectories
Run track_analysis_preprocess.m
Inputs: .csv file with all trajectories.
Outputs: .csv file with trajectories >= 10 frames and with a duty cycle >= 0.5.
Dependencies: none

4. Select and export trajectories of interest
Select trajectories for analysis within Swift. For SC-localized trajectories, use reference image of SC.
Calculate mean jump displacements and start times and duration if not previous analyzed.
Export as .csv
Inputs: filtered .csv file from 3 or raw trajectory information from 1 if not using prefilter.
Outputs: selected trajectories in .csv file
Dependencies: none

6. Run analysis software on curated trajectories.
Run track_analysis_main.m
Inputs: selected trajectories from 4.
Outputs: MSD plots, D values, and CDF of D values, including analysis for each principal component of motion.
Dependencies: processLags.m

