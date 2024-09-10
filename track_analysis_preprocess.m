%% 
% Lexy von Diezmann, 2023-2024. Released under the GNU GPL v3.
%
% This set of functions contains everything you need on the MATLAB side to process information 
% from Vutara SRX output files. It works in concert with the Swift program (Endesfelder 
% lab, using v0.43) to concatenate tracks.
%  
% Workflow:
% 

% See discussion e.g. in LvD and Rog JPCB 2021.
%
% This function

% Otherwise, start by analyzing raw localizations in Swift using standard settings (load each time). 
% Then, we will filter them here. Read in the csv:

currentTime = datetime('now','Format','yMMdd-HHmm');

[initialTracksName, initialTracksPath] = uigetfile('*.csv');

fullPath = [initialTracksPath initialTracksName];
T = readtable(fullPath);
numlocs=T.track_loc_count;
lifetime=T.track_lifetime;

figure; hist(numlocs./(lifetime+1),50)
valid = lifetime>=10;
valid = valid & numlocs./(lifetime+1)>=0.5;
T2=T;
T2(~valid,:) = [];

writetable(T2,[initialTracksPath char(currentTime) '_prefilter.csv'])
