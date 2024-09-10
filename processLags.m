% % Lexy von Diezmann, 2015-2024. Released under the GNU GPL v3.

% processLags: compute MSDs from arrays of displacements/lags
% 
% ALGORITHM SKETCH
% - Take displacements in large array packed with nans
% - find mean squared displacement as function of lag time, collapsing
% multiple dimensions into euclidean distance as needed
% - extract diffusion coefficient and dynamic error from linear fit using
% specified number of lags (could extend to anomalous diffusion or MLE in
% future but not currently supported)
% - output values of diffusion and error, and other intermediate values
%
% INPUTS
% "lagDisps" is an array of size #disps x #lags x #dimensions
% - #disps is the theoretical maximum number of displacements and is
% equal to the combined length of all tracks, but is packed with nans for
% i.e. cases where molecule blinked & as longer lags reduce number of total
% displacements with that lag measured in a given trajectory
% - #lags is the total number of lags calculated
% - #dimensions is the number of dimensions to take Euclidean norm of
% - to compute the average displacement for a given lag, take nanmean along
% each column. 
% - the units of lags should be in nm
% - In addition to cartesian coordinates, displacements 
% could be defined along cell axis (1D) or on the surface (2D)
%
% "lagTimes" is of size 1 x #lags and should be given in seconds
%
% "lagsToFit" defines the number of lags used to calculate MSD (typically
% 2-4: see e.g. Michalet and Berglund for discussion of tradeoffs)
%
% "frameTime" is exposure time in ms
%
% "model" defines the motion model to use. Current arguments:
% - "simpleLinear": 1D/2D/3D Brownian motion including motion blur
% - "surfaceLinear": as above but assume 3D motion is on 2D manifold
%
% "ID" is a binary array of size lagDisps that defines which lags to use
% (if e.g. selecting by pole, etc.)
%
% OUTPUTS
%
% "fitD": extracted diffusion coefficient in um^2/s
% "fitSigma": extracted dynamic error in nm
% "numDisps" number of displacements used to calculate MSD for each lag
% "MSD": mean square displacements for this set of lags
% "SDSD": standard deviation of square displacements 
% "MSDfit": structure containing fit values (slope, intercept)
%
%
function [fitD,fitSigma,numDisps,MSD,SDSD,MSDfit,bsD,bsSig,fitVals, fit1,fit2] = processLags(lagDisps,lagTimes,lagsToFit,frameTime,model,ID,resampNum,trackIDs)

if ~exist('model')
    model = 'simpleLinear';
end
if ~exist('ID')
    ID = true(size(lagDisps(:,:,1)));
end

expTime = frameTime;% / 1000; % exposure time in seconds
maxLags = size(lagDisps,2);
numDims = size(lagDisps,3);

lagEucDisps = sqrt(sum(lagDisps.^2,3)); % RMS disps (i.e. add in quadrature)
lagEucDisps(~ID) = nan; % strip out nonmembers as NaNs
%
numDisps = sum(~isnan(lagEucDisps),1);
%
SD = (lagEucDisps/1000).^2; % square displacements in um
MSD = nanmean(SD,1); % MSD in um^2, 1 x nLags
SDSD = nanstd(SD,[],1); % standard deviation of SDs, 1 x nLags

% fit the MSDs as a function of "lagTimes" for lags "lagsToFit"
% for linear fitting of isotropic brownian motion, 
% equations are MSD(tau,sigma) % tau: lag in seconds; sigma: dynamic errors
% MSD = 2nD*(tau-texp/3) + 2 sum (sigma^2)
% [ texp: exposure time in seconds ; n: #dimensions ; sum over dimensions]
% slope = MSD/tau = 2nD; D = slope / 2n, in um
% sigma = sqrt(intercept/2) + slope * texp / 6
% (not texp/3: see Backlund MBotC 2014 vs. Savin and Doyle 2005)
%
% for linear fits, structure out with .p1: slope. .p2: intercept.

MSDfit = fit(lagTimes(1:lagsToFit)',MSD(1:lagsToFit)','poly1');

if isequal(model,'surfaceLinear')&&numDims==3
    N = 2;
else
    N = numDims;
end

fit1 = MSDfit.p1;
fit2 = MSDfit.p2;
fitD = MSDfit.p1/(2*N);
fitSigma = sqrt(MSDfit.p2/(2*N)+fitD*expTime/3)*1000; % approximate sigma in nm

fitVals = MSDfit(lagTimes);

%% select specific tracks and set up bootstrapping
if exist('resampNum') && exist('trackIDs') && resampNum>1

% resample, keeping tracks pooled together
    
    % if it is desirable to save the index:
    % lagIdx = cell(resampNum,length(lagTimes));

trackIDs(~ID) = nan; % all track numbers in format of lagDisps;
uniqueTracks = unique(trackIDs(~isnan(trackIDs)));

bsD = nan(1,resampNum);

% bootstrap D values: define indexes for lags from track IDs
for bsIdx = 1:resampNum
    bsMSD = nan(size(MSD));
    for k = 1:maxLags
        tracksToUse = datasample(uniqueTracks,length(uniqueTracks));
%         tracksToUse = uniqueTracks; % this should yield same D as above
        idx = [];
        for l = 1:length(tracksToUse)
            idx = [idx;find(tracksToUse(l)==trackIDs(:,k))];
        end
%         lagIdx{bsIdx,k} = idx;
        bsMSD(k) = nanmean(SD(idx,k));        
    end
    bsFit = fit(lagTimes(1:lagsToFit)',bsMSD(1:lagsToFit)','poly1');    
    bsD(bsIdx) = bsFit.p1/(2*N);
    bsSig(bsIdx) = sqrt(bsFit.p2/(2*N)+bsD(bsIdx)*expTime/3)*1000;

    if bsIdx/100==round(bsIdx/100)
        disp(['reached iteration ' num2str(bsIdx)]);
    end
end 

else
    bsD = nan; % not computed
    bsSig = nan; % not computed
end % end resampling section
end % end function