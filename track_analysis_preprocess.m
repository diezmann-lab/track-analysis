%% 
% Lexy von Diezmann, 2023-2024. Released under the GNU GPL v3.

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
