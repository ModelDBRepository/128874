%% script to assess basic firing stats of simulations...

load('..\RandomInput.mat'); load('..\RandomInput_SIMPARAMS.mat')

% get spike data for neurons
MSspks = out.STms; 
MSspks(:,1) = MSspks(:,1)+1; % change from zero-base to 1-base index 
MSidxs = unique(MSspks(:,1));
Nmsidxs = numel(MSidxs);
Nms = SIMPARAMS.net.MS.N;

FSspks = out.STfs;
FSspks(:,1) = FSspks(:,1)+1; % change from zero-base to 1-base index 
FSidxs = unique(FSspks(:,1));
Nfsidxs = numel(FSidxs);
Nfs = SIMPARAMS.net.FS.N;

% simulation parameters
simT = SIMPARAMS.sim.tfinal; % in ms
T = simT * 1e-3;    % in seconds

% get stats for all MSNs
MSrate = zeros(Nms,1); MS_CVisi = zeros(Nms,1);
for j = 1:Nmsidxs
    currix = find(MSspks(:,1) == MSidxs(j));
    % basic rate
    MSrate(MSidxs(j)) = numel(currix) / T;
    
    % ISIs
    ts = MSspks(currix,2)*1e-3; % in seconds
    ISIs = diff(ts);
    MS_CVisi(MSidxs(j)) = std(ISIs)/mean(ISIs);
end
MS_CVisi(isnan(MS_CVisi)) = 0;

% get stats for all FSIs
FSrate = zeros(Nfs,1); FS_CVisi = zeros(Nfs,1);
for j = 1:Nfsidxs
    currix = find(FSspks(:,1) == FSidxs(j));
    % basic rate
    FSrate(FSidxs(j)) = numel(currix) / T;
    
    % ISIs
    ts = FSspks(currix,2)*1e-3; % in seconds
    ISIs = diff(ts);
    FS_CVisi(FSidxs(j)) = std(ISIs)/mean(ISIs);
end
FS_CVisi(isnan(FS_CVisi)) = 0;

%%%% plot empirical CDFs... 
%%% NOTE that stats structure returned by cdfplot has all info you need in
%%% it for all neurons
% with addition of:
MSrate_IQR = iqr(MSrate);
FSrate_IQR = iqr(FSrate);
MSCV_IQR = iqr(MS_CVisi);
FSCV_IQR = iqr(FS_CVisi);

figure(10); clf
[hM,MSstats] = cdfplot(MSrate); hold on
[hF,FSstats] = cdfplot(FSrate); set(hF,'Color',[1 0 0])
phF = get(hF,'Parent'); set(phF,'XScale','log')
title('ECDF plots for firing rates')
xlabel('x=firing rate (spikes/s)')

figure(11); clf
[hMCV,MSCVstats] = cdfplot(MS_CVisi); hold on
[hFCV,FSCVstats] = cdfplot(FS_CVisi(~isnan(FS_CVisi))); set(hFCV,'Color',[1 0 0])
title('ECDF plots for ISI CV')
xlabel('x=ISI CV')

n_active_MSNs = sum(MSrate > 0)
n_active_FSIs = sum(FSrate > 0)


