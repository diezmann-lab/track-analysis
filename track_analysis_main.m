%% 
% Lexy von Diezmann, 2023-2024. Released under the GNU GPL v3.

fileSuffix = '_final'; % adjust manually
currentTime = datetime('now','Format','yMMdd-HHmm');
fileSuffix = [fileSuffix char(currentTime)];
[selectedTracksName, selectedTracksPath] = uigetfile('*.csv','MultiSelect','on');

if ~iscell(selectedTracksName)
    selectedTracksName = {selectedTracksName};
end
dataAll = cell(size(selectedTracksName));
tracksAll = cell(size(dataAll));

% try to include all these variables; remove any that are not present
testNames = {'x','y','z','t','id','precisionx','precisiony','precisionz','seg_id','track_id'};
for i = 1:length(selectedTracksName) 
    fullPath = [selectedTracksPath selectedTracksName{i}];
    opts = detectImportOptions(fullPath);
    
    allNames = opts.SelectedVariableNames;
    testNames = intersect(testNames,allNames);
    opts.SelectedVariableNames = testNames;
    clear allNames
    
    dataAll{i} = readtable(fullPath,opts);
    
end
clear testNames dataIn

%% analyze all files

for fileNum = 1:length(dataAll)
    data = dataAll{fileNum};

%% More preprocessing and tabulation

frameTime = 0.1; % seconds
minLocNum = 10; % minimum number to include in analysis

trackcounter = 1;
segcounter = 1;

% swift outputs both "track" and sub-track "segment" annotations
segid = data.seg_id;
trackid = data.track_id;

segment = nan(size(segid));
track = nan(size(trackid));
segsToUse = unique(segid)';
tracksToUse = unique(trackid)';
% segsToUse(ismember(segsToUse,fastsegs)) = [];

for i = segsToUse
    segment(segid==i)=segcounter;
    segcounter=segcounter+1;
end

for i = tracksToUse
    track(trackid==i)=trackcounter;
    trackcounter=trackcounter+1;
end

% filter based on track or segment? Adjust manually.
options = 'track';
clear trackid segid segsToUse tracksToUs

%% process tracks, including making continuous by filling in NaNs

x = data.x;
y = data.y;
t = data.t;
id = data.id; % unique identifier for each localization

trackInfo = struct;

nLags = 40;
lagTimes = (1:nLags)*frameTime;
lagsToFit = 4;
msdModel = 'simpleLinear';
useBS = 0;


switch options
    case 'track'
        splitBy = track;
    case 'segment'
        splitBy = segment;
    otherwise
end

% for loop populating track by track (or segment by segment)
for i = 1:max(splitBy)
    xLoc = x(splitBy==i);
    yLoc = y(splitBy==i);
    frame = t(splitBy==i);
    molID = id(splitBy==i);
    xDiffs = []; yDiffs = []; frameDiffs = [];
    while length(frame)<(frame(end)-frame(1)+1)
        frameDiff = diff(frame)-1;
        idx = find(frameDiff>0,1);
        xLoc = [xLoc(1:idx); nan(frameDiff(idx),1); xLoc(idx+1:end)];
        yLoc = [yLoc(1:idx); nan(frameDiff(idx),1); yLoc(idx+1:end)];
        % zLoc = [zLoc(1:idx); nan(frameDiff(idx),1); zLoc(idx+1:end)];
        frame = [frame(1:idx); nan(frameDiff(idx),1); frame(idx+1:end)];
        molID = [molID(1:idx); nan(frameDiff(idx),1); molID(idx+1:end)];
    end
    for k = 1:min(nLags,length(frame)-1)
        xDiffs(:,k) = [xLoc(1+k:end)-xLoc(1:end-k) ; nan(k-1,1)];
        yDiffs(:,k) = [yLoc(1+k:end)-yLoc(1:end-k) ; nan(k-1,1)];
        frameDiffs(:,k) = [frame(1+k:end)-frame(1:end-k) ; nan(k-1,1)];
    end
    if k<nLags
        xDiffs(:,k+1:nLags) = nan;
        yDiffs(:,k+1:nLags) = nan;
        frameDiffs(:,k+1:nLags) = nan;
    end

    mjd = nanmean(sqrt(xDiffs(:,1).^2+yDiffs(:,1).^2)); 
    msd = nanmean(xDiffs.^2+yDiffs.^2)/(1000^2); % convert to um. 
    % note that this is not the same: mjd = sqrt(msd(:,1));
  
    % generate per-track MSD. note this can fail with short tracks

    if any(isnan(msd(1:lagsToFit)))
        fitD = nan;
        fitSigma = nan;
    else
        MSDfit = fit(lagTimes(1:lagsToFit)',msd(1:lagsToFit)','poly1');
        N = 2; %2-dimensional MSD
        fit1 = MSDfit.p1;
        fit2 = MSDfit.p2;
        fitD = MSDfit.p1/(2*N);
        fitSigma = sqrt(MSDfit.p2/(2*N)+fitD*frameTime/3)*1000; % approximate sigma in nm
    end
    trackInfo(i).xLoc = xLoc;
    trackInfo(i).yLoc = yLoc;
    trackInfo(i).frame = frame;
    trackInfo(i).molID = molID;
    trackInfo(i).xDiffs = xDiffs;
    trackInfo(i).yDiffs = yDiffs;
    trackInfo(i).frameDiffs = frameDiffs;
    trackInfo(i).mjd = mjd;
    trackInfo(i).msd = msd;
    trackInfo(i).length = length(xLoc);
    trackInfo(i).numLocs = sum(~isnan(xLoc));
    trackInfo(i).fitD = fitD;
    trackInfo(i).fitSigma = fitSigma;
end

tracksAll{fileNum} = trackInfo;

end

trackInfo = [tracksAll{:}];

%% do PCA on tracks including MSD along principal components

for i = 1:length(trackInfo)
    xy = [trackInfo(i).xLoc-mean(trackInfo(i).xLoc,"omitnan"),...
          trackInfo(i).yLoc-mean(trackInfo(i).yLoc,"omitnan")];    
    pMat = pca(xy);
    horzXY = xy*pMat;

    p1diffs = trackInfo(i).xDiffs*pMat(1,1) + trackInfo(i).yDiffs*pMat(2,1);
    p2diffs = trackInfo(i).xDiffs*pMat(1,2) + trackInfo(i).yDiffs*pMat(2,2);

    msd1 = nanmean(p1diffs.^2)/(1000^2); % convert to um
    msd2 = nanmean(p2diffs.^2)/(1000^2); % convert to um
  
    % generate per-track MSD. note this often fails with jumps/short tracks
    
    % primary dimension
    lagsToFit = 4;
    if any(isnan(msd1(1:lagsToFit)))
        fitD1 = nan;
        fitSigma1 = nan;
    else
        MSDfit = fit(lagTimes(1:lagsToFit)',msd1(1:lagsToFit)','poly1');
        N = 1; %1-dimensional MSD
        fit1 = MSDfit.p1;
        fit2 = MSDfit.p2;
        fitD1 = MSDfit.p1/(2*N);
        fitSigma1 = sqrt(MSDfit.p2/(2*N)+fitD1*frameTime/3)*1000; % approximate sigma in nm
    end
    % secondary dimension
    lagsToFit = 4;
    if any(isnan(msd2(1:lagsToFit)))
        fitD2 = nan;
        fitSigma2 = nan;
    else
        MSDfit = fit(lagTimes(1:lagsToFit)',msd2(1:lagsToFit)','poly1');
        N = 1; %1-dimensional MSD
        fit1 = MSDfit.p1;
        fit2 = MSDfit.p2;
        fitD2 = MSDfit.p1/(2*N);
        fitSigma2 = sqrt(MSDfit.p2/(2*N)+fitD2*frameTime/3)*1000; % approximate sigma in nm
    end
    dim1std = std(horzXY(:,1),'omitnan');
    dim2std = std(horzXY(:,2),'omitnan');
    anisotropy = dim1std/dim2std;
    trackInfo(i).dim1std = dim1std;
    trackInfo(i).dim2std = dim2std;
    trackInfo(i).anisotropy = anisotropy;
    trackInfo(i).pMat = pMat;
    trackInfo(i).p1 = horzXY(:,1);
    trackInfo(i).p2 = horzXY(:,2);
    trackInfo(i).p1diffs = p1diffs;
    trackInfo(i).p2diffs = p2diffs;
    trackInfo(i).fitD1 = fitD1;
    trackInfo(i).fitD2 = fitD2;
    % trackInfo(i).fitSigma1 = fitSigma1; 
    % trackInfo(i).fitSigma2 = fitSigma2; % often complex
    clear pMat xy horzXY dim1std dim2std anisotropy p1diffs p2diffs
end

%% 
% Let's try looking at those D values on a per-track basis.
%

minLocNum = 10;

allDs = [trackInfo(:).fitD];
allNumLocs = [trackInfo(:).numLocs];

goodTracks = allNumLocs>minLocNum & allDs>0;

weightedD = [];
weightedMJD = [];
if isfield(trackInfo,'fitD1')
    weightedD1 = [];
    weightedD2 = [];
end

totalLocs = 0;
for i = find(goodTracks)
    weightedD = [weightedD repmat(trackInfo(i).fitD,1,allNumLocs(i))];
    if isfield(trackInfo,'fitD1')
        weightedD1 = [weightedD1 repmat(trackInfo(i).fitD1,1,allNumLocs(i))];
        weightedD2 = [weightedD2 repmat(trackInfo(i).fitD2,1,allNumLocs(i))];
    end
    totalLocs = totalLocs + allNumLocs(i);
end

% figure; hist([trackInfo(goodTracks).fitD],30);
binSize = 0.0025;

[vals,x] = hist(weightedD(weightedD>0),binSize/2:binSize:0.05);
vals=vals/sum(vals);
h = figure;
% bar(x,vals);
plot(x,vals);

xlabel('estimated D value (um2/s)'); ylabel('number of localizations in tracks with this value')
xlim([0 0.05]);

saveas(h,['histogram_plot' fileSuffix],'fig'); close(h);

% output cumulative distribution functions

h = figure; ecdf(weightedD(weightedD>0));
xlim([0 0.1]); title(['D cdf ' fileSuffix]);
saveas(h,['D_ecdf' fileSuffix],'fig'); close(h);

h = figure; ecdf(weightedD1(weightedD1>0));
xlim([0 0.1]);title(['D1 cdf ' fileSuffix]);
saveas(h,['D1_ecdf' fileSuffix],'fig'); close(h);

h = figure; ecdf(weightedD2(weightedD2>0));
xlim([0 0.1]);title(['D2 cdf ' fileSuffix]);
saveas(h,['D2_ecdf' fileSuffix],'fig'); close(h);
clear vals x


%% plot all selected tracks

numLocs = [trackInfo(:).numLocs];
trackLength = [trackInfo(:).length];
allDs = [trackInfo(:).fitD];
allDs1 = [trackInfo(:).fitD1];

% cutoff=50;

idx = find(numLocs>70 & allDs>0 & allDs1<0.03); % setting for ZHP-3 EP
idx = find(numLocs>52 & allDs>0 & allDs1<0.03); % setting for ZHP-3 LP
idx = find(numLocs>70 & allDs>0 & allDs1<0.03); % setting for SYP-3 EP
% idx = find(numLocs>110 & allDs>0 & allDs1<0.03); % setting for SYP-3 LP
% idx = find(numLocs>10 & allDs>0);
% idx = find(trackLength>70 & allDs>0 & allDs<0.02);
tracksToUse = trackInfo(idx);


% do a sort
dVals = [tracksToUse.fitD];

dVals = [tracksToUse.fitD1];

numLocs = [tracksToUse.numLocs];
trackLength = [tracksToUse.length];
anisotropy = [tracksToUse.anisotropy];
dim1sigma = [tracksToUse.dim1std];
dim2sigma = [tracksToUse.dim2std];
% [Dvals,sortedIdx] = sort(Dvals);

[~,sortedIdx] = sort(anisotropy);
[~,sortedIdx] = sort(dVals);
% [~,sortedIdx] = sort(dim1sigma);
tracksToUse = tracksToUse(sortedIdx);

% regenerate values for sorted tracks
dVals = [tracksToUse.fitD];


% for D1
dVals = [tracksToUse.fitD1];

numLocs = [tracksToUse.numLocs];
trackLength = [tracksToUse.length];
anisotropy = [tracksToUse.anisotropy];
dim1sigma = [tracksToUse.dim1std];
dim2sigma = [tracksToUse.dim2std];
%
rangeMax = max(dVals);
rangeMax = 0.03;

% can change colors to indicate other parameters
% colors = colormap(parula(length(sortedIdx)));

numColors = 256;

dRange = linspace(0,rangeMax,numColors);

colorIdx = round(interp1(dRange,1:numColors,dVals));

colors = parula(numColors);

censored = []; % can mark some tracks as gray to highlight others
% censored = find(dVals<0.002 | dVals>0.02);
h = figure; hold on;

allX = [];
allY = [];

for i = 1:length(tracksToUse)

x = tracksToUse(i).p2; x(isnan(x)) = [];
y = tracksToUse(i).p1; y(isnan(y)) = [];

if ~ismember(i,censored)
    allX = [allX;x];
    allY = [allY;y];
end
shiftRow = 32; % number of tracks before shifting
xOffset = 500; % between tracks in nm
yOffset = 2400; % between rows in nm
% xOffset = 0; yOffset = 0; % to plot all tracks in same coord plane
% by default, MATLAB uses a high-contrast repeating set of 7 colors
% plot(x+mod(i,shiftRow)*xOffset,y-yOffset*floor(i/shiftRow),'.-')

% if you specified "colors" above, could use that instead

% can censor as gray if defined
    if ismember(i,censored)
        plotColor = [0.5 0.5 0.5];
    elseif exist('colorIdx')
        plotColor = colors(colorIdx(i),:);
    end
    plot(x+mod(i,shiftRow)*xOffset,y-yOffset*floor(i/shiftRow),'.-','color',plotColor)

end
axis equal

cbh = colorbar;

clim([dRange(1),dRange(end)]);
cbh.Limits = [dRange(1),dRange(end)];

title('example tracks');

if exist('fileSuffix')
    saveas(h,['aligned tracks' fileSuffix],'fig');
    saveas(h,['aligned tracks' fileSuffix],'epsc'); close(h);
end

%% Example code to plot "averaged track" in shared coordinate system.
% 
% % plot X profile
% 
% figure; plot(allX,allY,'.');
% 
%     range = linspace(-250,250,100);
%     vals = hist(allX,range);
%     figure; bar(range,vals); xlim([-250 250]);

%% do averaged MSD prep

allLagsX = nan(0,nLags);
allLagsY = nan(0,nLags);
allLagIDs = nan(0,nLags);
allLengths = [trackInfo(:).length];

numLocs = [trackInfo(:).numLocs];
allDs = [trackInfo(:).fitD];
allD1s = [trackInfo(:).fitD1];

% If the prefilter step is skipped, tracks can be filtered for length and duty cycle here.
idx = find(numLocs>10 & allDs>0);
tracksToUse = trackInfo(idx);

medianLength = median([tracksToUse(:).length]);
medianDuty = median([tracksToUse(:).numLocs]./[tracksToUse(:).length]);
for i = idx
    allLagsX = [allLagsX; trackInfo(i).xDiffs];
    allLagsY = [allLagsY; trackInfo(i).yDiffs];
    allLagIDs = [allLagIDs; i*ones(size(trackInfo(i).xDiffs))];
end
MSDx = nanmean(allLagsX.^2,1);
MSDy = nanmean(allLagsY.^2,1);

if isfield(trackInfo,'p1diffs')
    allLags1 = nan(0,nLags);
    allLags2 = nan(0,nLags);
    for i = idx
        allLags1 = [allLags1; trackInfo(i).p1diffs];
        allLags2 = [allLags2; trackInfo(i).p2diffs];
    end
end
MSDxy = MSDx+MSDy;

lagTimes = (1:nLags)*frameTime;

%% do MSD calc for averaged data

lagsToFit = 4;
msdModel = 'simpleLinear';
allLagsXY = cat(3,allLagsX,allLagsY);

% use bootstrapping to estimate errors? (Manually adjust)
useBS = 1; bsNum = 50; % number of realizations
% useBS = 0;

numTracksAnalyzed = length(unique(allLagIDs(:)));

if useBS
    [allD,allSig,allNumDisps,allMSD,allSDSD,~,allDbs,allSigbs] = ...
        processLags(allLagsXY,lagTimes,lagsToFit,frameTime,msdModel,ones(size(allLagsXY)),bsNum,allLagIDs);
    h = figure; errorbar(lagTimes,allMSD,allSDSD./sqrt(allNumDisps));
    xlim([0 2]);
    xlabel('lag time (s)'); ylabel('MSD (um^2)');
    title(['MSD plot: D = ' num2str(allD*1000,4) ' \pm ' num2str(std(allDbs*1000),4) ' x 10^{-3 um^2/s}, ', num2str(numTracksAnalyzed) ' unique tracks']);
else
    [allD,allSig,allNumDisps,allMSD,allSDSD,~,~,~,fitVals] = ...
        processLags(allLagsXY,lagTimes,lagsToFit,frameTime,msdModel,ones(size(allLagsXY)));
    h = figure; errorbar(lagTimes,allMSD,allSDSD./sqrt(allNumDisps));
    xlim([0 2]);
    xlabel('lag time (s)'); ylabel('MSD (um^2)');
    title(['MSD plot: D = ' num2str(allD*1000,4) ' x 10^{-3} µm^2/s']);
end
saveas(h,['MSD_curve' fileSuffix],'fig');
saveas(h,['MSD_curve' fileSuffix],'epsc');
close(h);


if isfield(trackInfo,'p1diffs')
    % P1 MSD
    if useBS % P1 MSD with bootstrapping
    [allD1,allSig1,allNumDisps,allMSD,allSDSD,~,allDbs1,allSigbs] = ...
        processLags(allLags1,lagTimes,lagsToFit,frameTime,msdModel,ones(size(allLags1)),bsNum,allLagIDs);
    h = figure; errorbar(lagTimes,allMSD,allSDSD./sqrt(allNumDisps));
    xlim([0 2]);
    xlabel('lag time (s)'); ylabel('MSD (um^2)');
    title(['MSD plot: D = ' num2str(allD1*1000,4) ' \pm ' num2str(std(allDbs1*1000),4) ' x 10^-3 um^2/s']);
    else % P1 MSD without bootstrapping
    [allD1,allSig1,allNumDisps,allMSD,allSDSD] = ...
        processLags(allLags1,lagTimes,lagsToFit,frameTime,msdModel,ones(size(allLags1)));
    h = figure; errorbar(lagTimes,allMSD,allSDSD./sqrt(allNumDisps));
    xlim([0 2]);
    xlabel('lag time (s)'); ylabel('MSD (um^2)');
    title(['MSD plot: D (P1) = ' num2str(allD*1000,4) ' x 10^{-3} µm^2/s']);
    end
    saveas(h,['MSD_curve_P1' fileSuffix],'fig');
    saveas(h,['MSD_curve_P1' fileSuffix],'epsc'); close(h);

    if useBS % P2 MSD with bootstrapping
    [allD2,allSig2,allNumDisps,allMSD,allSDSD2,~,allDbs2,allSigbs] = ...
        processLags(allLags2,lagTimes,lagsToFit,frameTime,msdModel,ones(size(allLags2)),bsNum,allLagIDs);
    h = figure; errorbar(lagTimes,allMSD,allSDSD./sqrt(allNumDisps));
    xlim([0 2]);
    xlabel('lag time (s)'); ylabel('MSD (um^2)');
    title(['MSD plot: D = ' num2str(allD2*1000,4) ' \pm ' num2str(std(allDbs2*1000),4) ' x 10^-3 um^2/s']);
    else % P2 MSD without bootstrapping

    [allD2,allSig2,allNumDisps,allMSD,allSDSD,~,allDbs2,allSigbs] = ...
        processLags(allLags2,lagTimes,lagsToFit,frameTime,msdModel,ones(size(allLags2)));
    h = figure; errorbar(lagTimes,allMSD,allSDSD./sqrt(allNumDisps));
    xlim([0 2]);
    xlabel('lag time (s)'); ylabel('MSD (um^2)');
    title(['MSD plot: D (P2) = ' num2str(allD2*1000,4) ' x 10^{-3} µm^2/s']);
    end
    saveas(h,['MSD_curve_P2' fileSuffix],'fig');
    saveas(h,['MSD_curve_P2' fileSuffix],'epsc'); close(h);
end

save(['MSD_analyze' fileSuffix '.mat'],'allD','allSig','allMSD','allSDSD','allDbs','allSigbs','trackInfo','tracksToUse','allLagsX','allLagsY', 'allLags1','allLags2','lagTimes','lagsToFit','frameTime', 'allLagIDs', 'allNumDisps','frameTime','weightedD','allD1','allDbs1','allD2','allDbs2','numTracksAnalyzed','medianLength','medianDuty')
