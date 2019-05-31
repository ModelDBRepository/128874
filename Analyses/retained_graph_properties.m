%%% post-hoc analysis of best clustered binsize, run after
%%% cluster_words_winsize
clear all
load('first_stage_analysis')


%---------------- firing stats of whole network --------------------------
firing_stats

%------------------ effect of changing analysis parameters ---------------
% can check effect of changing thetaH on the best one!!
% MAKE SURE YOU REGENERATE hc BEFORE DOING THIS....
HE = pdist(scount{btest}','hamming');
hc = squareform(HE);    % turn it into matrix - with zeroes on diagonal

thetaHvec = 0.05:0.05:0.3;
rtnd = cell(numel(thetaHvec),1);
Sc_thetaH = cell(numel(thetaHvec),1);
CmplxH = zeros(numel(thetaHvec),1);

for i = 1:numel(thetaHvec)
    i
    Ac = (hc <= thetaHvec(i));    % smallest error
    dgrs = sum(Ac);
    delN = find(dgrs < 2);
    rtnd{i} = setdiff(idxs,delN);    %% index back into Rn to recover original index of MSN...

    % get reduced version (all nodes degree at least 2 before reduction...)
    Acs = Ac(rtnd{i},rtnd{i});
    edgs(i) = sum(sum(Acs)); nds(i) = numel(rtnd{i});
    if nds(i) > nlimit & edgs(i) > log(nds(i))   % giant component phase transition 
        % only group if there are enough nodes and  edges to make it worthwhile...
        
        [Sc_thetaH{i},Uc] = multileadevsplit(Acs) ;% ,'refine'); % do refine step later when got time...
        thetaHngrps(i) = max(Sc_thetaH{i});
        siz = [];
        for j = 1:thetaHngrps(i)
            siz = [siz numel(find(Sc_thetaH{i} == j))];   
        end
        thetaHgrpsizes{i} = siz;
    else
        Sc_thetaH{i} = [];  % nothing to group!
        thetaHngrps(i) = 0;
        thetaHgrpsizes{i} = [];
    end
    %%% Hamming distance within groups vs whole retained network for
    %%% modularity method
    Hamming_all_H(i) = median(pdist(scount{btest}(:,rtnd{i})','hamming'));
    for j=1:thetaHngrps(i)
        cts = [scount{btest}(:,rtnd{i}(Sc_thetaH{i}==j))];
        Hamming_H{i}(j) = median(pdist(cts','hamming'));
    end

    % note that complexity measure now has same Hdiff all the time, so
    % changes only with changes in numberof groups and proportion of neuron
    if thetaHngrps(i)
        CmplxH(i) = sum(thetaHgrpsizes{i})/Nms * thetaHngrps(i)  * hdiff(btest);
    end
end

% and threshold group size of interest
grpthresh = round(Nms * 0.01);
thetaHngrpsT = zeros(numel(thetaHngrps),1);
for i = 1:numel(thetaHngrps)
    thetaHngrpsT(i) = sum(thetaHgrpsizes{i} >= grpthresh);
end

%%%% plot out some results...
figure(99), clf
subplot(211), plot(thetaHvec,thetaHngrps,'+-'), ylabel('# grps'), xlabel('Threshold')
subplot(212), plot(thetaHvec,CmplxH,'+-'), ylabel('\beta'),  xlabel('Threshold')
%plot(thetaHvec,thetaHngrpsT,'k+-')

% find most groups
% bH = find(thetaHngrps == max(thetaHngrps)); 
% bH = find(thetaHngrpsT == max(thetaHngrpsT));
bH = find(CmplxH == max(CmplxH));
if numel(bH) > 1 bH = bH(1); end
%bH = 2;

ScH = Sc_thetaH{bH};
RH = rtnd{bH};
Ct = scount{btest};

number_of_groups = thetaHngrps(bH)

grp_membersH = {};
for loop = 1:thetaHngrps(bH)
   grp_membersH{loop} = RH(ScH==loop); 
end

bins = 0:binsize(btest):simT;
sctsH = zeros(numel(bins),thetaHngrps(bH));
figure(103), clf
for loop = 1:thetaHngrps(bH)
    cts = [Ct(:,RH(ScH==loop))]';
    sctsH(:,loop) = [sum(cts)./thetaHgrpsizes{bH}(loop)]';
    subplot(max(ScH),1,loop), bar(bins,sctsH(:,loop),'histc')
    % subplot(max(ScH),1,loop), bar(smcts(:,loop))
    title('Summed spike count for graph-method')
end

figure(104), clf, hold on, %title('Grouped raster by graph-method (group 1 at the top)')
% full colour version
clrs = colormap; [rw cl] = size(clrs); clrstep = floor(rw / thetaHngrps(bH));
% alternating black and grey version
%clrs = repmat([0 0 0; 0.5 0.5 0.5],thetaHngrps(bH),1); clrstep = 1;

ctr = numel(RH);
for loop = 1:thetaHngrps(bH)
    inds = RH(ScH==loop);
    %inds = RH(ScH == grp_seq_hand(loop));
    stgrp = [];
    ixgrp = [];
    for i = 1:numel(inds)
        indst = find(MSspks(:,1) == inds(i));
        st = MSspks(indst,:);
        % plot(st(:,2), ones(length(st),1)+ctr,'.','Color',clrs(clrstep*loop,:));
        % plot(st(:,2), ones(length(st),1)+ctr,'.');
        
        stgrp = [stgrp; st(:,2)];
        ixgrp = [ixgrp; ones(numel(st(:,1)),1)+ctr];
        ctr = ctr-1;
    end
    plot(stgrp, ixgrp,'.','Color',clrs(clrstep*loop,:));
    %fname1 = ['ph_stgrp' num2str(loop) '.txt'];  fname2 = ['ph_ixgrp' num2str(loop) '.txt']; 
    %save(fname1,'stgrp','-ascii'); save(fname2,'ixgrp','-ascii');
end
% axis off

%%% plot unsorted raster to drive home how good the clustering is...
figure(105), clf, hold on, %title('Ungrouped raster of all clustered cells')
for loop = 1:numel(rtnd{bH})
    ixs = find(MSspks(:,1) == rtnd{bH}(loop));
    plot(MSspks(ixs,2),ones(numel(ixs),1)*loop,'k.');
end
% axis off

%%%% now go through and get index into best stuff from this analysis

%------------------ firing properties of each cluster --------------------
%%% get firing rates and spectra of spike trains in clusters....
grprates = cell(thetaHngrps(bH),1);
for loop = 1:thetaHngrps(bH)
   inds = grp_membersH{loop};
   for i = 1:thetaHgrpsizes{bH}(loop)
      currix = find(MSspks(:,1) == MSidxs(i));
      ts = MSspks(currix,2)*1e-3; % in seconds
      % difficult to compute spectrum cos of so few spikes in each cell...
   end
   grprates{loop} = MSrate(inds);
   grpmedian(loop) = median(MSrate(inds));
end



% find clusters that are antagonisitc
cxy = corrcoef(sctsH); % compute correlation coefficient between total spike counts...

% firing sequences?? Here looking at sequences of *silence* 
% 1. by min active vector in bin
TsctsH = sctsH(1:end-1,:);
min_in_bin = min(TsctsH');    % find min proportion active in each bin
grp_seq = zeros(numel(min_in_bin),1);
for i = 1:numel(min_in_bin)
    [x,grp_seq(i)] = ind2sub(size(TsctsH),find(TsctsH == min_in_bin(i),1));
end

% 2. by all vectors below some threshold in bin
grp_seq2 = cell(numel(min_in_bin),1);
for i = 1:numel(min_in_bin)
    grp_seq2{i} = find(TsctsH(i,:) < 0.2); 
end

%---------------- reconstruct correlation graph of best theta
Ac = (hc <= thetaHvec(bH));    % smallest error

% get reduced version (all nodes degree at least 2 before reduction...)
Acs = Ac(rtnd{bH},rtnd{bH});
Kdist = sum(Acs);

%----------------------------- properties of correlation graph  ---------------------------

% a few graph properties
Kmn = mean(Kdist); nds = numel(Kdist); dnsty = sum(Kdist) / (nds*(nds-1));
LR = log(nds) / log(Kmn);
CR = Kmn / nds;

%%----------------------------- do graph analysis of network --------------
%%%% get path-lengths of just MS-MS network using Olaf's code (reachdist.m)
Cmsms = SIMPARAMS.net.Cmsms +1; % was zero-indexed
Cmsms_b = SIMPARAMS.net.Cmsms_b+1; % was zero-indexed
Amsms = zeros(Nms);

% construct adjacency matrix for MS-MS connections
for i = 1:Nms
    % loop round and stick each MSN's connections into matrix....
    cur = Cmsms([Cmsms_b(i):Cmsms_b(i+1)-1]);
    Amsms(i,cur) = 1;
end

% properties
indgr = sum(Amsms);
outdrg = sum(Amsms');
Mmsms = sum(sum(indgr)); mK = Mmsms/Nms;  % approx mean degree...
msms_dnsty = Mmsms / (Nms*(Nms-1));

% get modularity!!!
[Smsms,U] = multileadevsplit(Amsms);

% get path lengths
[R,D] = reachdist(Amsms);
L = sum(sum(D))/Nms^2;    LR = log(Nms) / log(mK);
Ctri = clusttriang(Amsms);  CR = mK / Nms;
[msmsCws,Cws] = clustind(Amsms);
St = (Ctri/CR) / (L/LR);
Sws = (Cws/CR) / (L/LR); 

%--- go through groups and get within group MSN-MSN path length, clustering coeff..
grpL = zeros(thetaHngrps(bH),1);
grpCt = zeros(thetaHngrps(bH),1);
grpMC_WS = zeros(thetaHngrps(bH),1);
grpC_WS = cell(thetaHngrps(bH),1);
grpSCtri = zeros(thetaHngrps(bH),1);
grpSCws = zeros(thetaHngrps(bH),1);

figure(101), clf, hold on
for i = 1:thetaHngrps(bH)
    cur = grp_membersH{i};
    nmem = thetaHgrpsizes{bH}(i);
    curdist = D(cur,cur);
    % grpL(i) = median(curdist(:));
    % path lengths and clustering coefficients
    grpL(i) = sum(sum(D(cur,cur)))/nmem^2;
    grpCt(i) = clusttriang(Amsms(cur,cur));
    [grpC_WS{i},grpMC_WS(i)] = clustind(Amsms(cur,cur));
    % small-world-ness
    grp_m = sum(sum(Amsms(cur,cur)));   % number of edges
    mdeg =grp_m/ nmem;    % i.e. mean number of edges per node 
    % approximate random graphs values
    grpCRa = mdeg / nmem; grpCRb = mdeg / nmem;    % <k> / n
    grpLR = log(nmem) / log(mdeg); % log(n) / log <k>
    % Monte Carlo random graph values... [these groups are very small!!
%     grpLR = 0; grpCRa = 0; grpCRb = 0; nMC = 20;
%     for j = 1:nMC
%         % create random_graph of same number of edges
%         Arand = random_graph(nmem,0,grp_m,'directed');
%         % compute path length
%         [temp,tempD] = reachdist(Arand);
%         grpLR = grpLR + (sum(sum(tempD)) / nmem^2)/nMC;
%         % compute mean clustering coeffs
%         grpCRa = grpCRa + clusttriang(Arand)/nMC;
%         [temp,tempC] = clustind(Arand);
%         grpCRb = grpCRb + tempC/nMC;
%     end
    
    grpSCtri(i) = (grpCt(i)/grpCRa) / (grpL(i)/grpLR);
    grpSCws(i) = (grpMC_WS(i)/grpCRb) / (grpL(i)/grpLR);
    % plot ECDF of clustering coefficient distribution...
    subplot(211), hold on, if ~isempty(grpC_WS{i}) [h,stats] = cdfplot(grpC_WS{i}); end
    set(h,'Color',clrs(clrstep*i,:)); xlabel('x = clustering coefficient')
    % plot histogram of path length distribution in each group
    subplot(212), hold on, plot(hist(curdist(:),max(curdist(:)))./nmem.^2,'+-','Color',clrs(clrstep*i,:));
    xlabel('Path length within group'); ylabel('probability of path length')
end
Hgrps = Hamming{bH}';


%--- go through all anatagonistic group pairs and get path length, clustering
%etc...
c_up = triu(cxy);   % just upper triangular part as correlation matrix is symmetric
ant_pairs = find(c_up <= -0.1);
[r,c] = ind2sub(size(c_up),ant_pairs);
antL_AB = zeros(numel(ant_pairs),1);
antL_BA = zeros(numel(ant_pairs),1);

% colormap = jet;
% clrs = colormap; [rw cl] = size(clrs); antclrstep = floor(rw / numel(ant_pairs));
figure(102), clf, hold on
for i = 1:numel(ant_pairs)
    curA = grp_membersH{r(i)};
    curB = grp_membersH{c(i)};
    nmem = numel(curA)*numel(curB);
    curdistAB = D(curA,curB);
    curdistBA = D(curB,curA);
    % antL(i) = median(curdist(:));
    %antL(i) = sum(sum(D(curR,curC)))/nmem;
    antL_AB(i) = mean(curdistAB(:));
    antL_BA(i) = mean(curdistBA(:));
    plot(hist(curdist(:),max(curdist(:)))./nmem,'+-','Color',clrs(antclrstep*i,:)); %% something broken here??
    xlabel('Path length between groups'); ylabel('probability of path length')
end


%----- construct adjacency matrix for FS-MS connections
Nfs = SIMPARAMS.net.FS.N;
Cfsms = SIMPARAMS.net.Cfsms +1; % was zero-indexed
Cfsms_b = SIMPARAMS.net.Cfsms_b+1; % was zero-indexed

Afsms = zeros(Nfs,Nms);
for i = 1:Nfs
    % loop round and stick each FSN's connections into matrix....
    cur = Cfsms([Cfsms_b(i):Cfsms_b(i+1)-1]);
    Afsms(i,cur) = 1;
end

% check that all MSNs get at least one FS input..
fsmsin = sum(Afsms);

% get FS indexes for each group...
grpFS = zeros(Nfs,thetaHngrps(bH));
grpFSw = zeros(Nfs,thetaHngrps(bH));
for i = 1:thetaHngrps(bH)
    cur = grp_membersH{i};
    grpFS(:,i) = [sum(Afsms(:,cur)')]'./thetaHgrpsizes{bH}(i);
    % weight by firing rate of FS interneuron
    grpFSw(:,i) = (FSrate .* [sum(Afsms(:,cur)')]')./thetaHgrpsizes{bH}(i);
end

% compare this to cxy - correlation coefficient between spike counts?

% construct adjacency matrix for FS-FS connections
Cfsfs = SIMPARAMS.net.Cfsfs +1; % was zero-indexed
Cfsfs_b = SIMPARAMS.net.Cfsfs_b+1; % was zero-indexed
Afsfs = zeros(Nfs);
for i = 1:Nfs
    % loop round and stick each MSN's connections into matrix....
    cur = Cfsfs([Cfsfs_b(i):Cfsfs_b(i+1)-1]);
    Afsfs(i,cur) = 1;
end

%%% Cgapfs indexes into Pgap matrix of connected pairs... 
Pgapfs = SIMPARAMS.net.Pgapfs +1; % was zero-indexed
Agapfs = zeros(Nfs);
for i = 1:length(Pgapfs)
    % loop round and stick each MSN's connections into matrix....
    Agapfs(Pgapfs(i,1),Pgapfs(i,2)) = 1;
    Agapfs(Pgapfs(i,2),Pgapfs(i,1)) = 1;
end

%%% whole network??
Astr = zeros(Nfs+Nms);
Astr(Nfs+1:end,Nfs+1:end) = Amsms;  % add MS MS connections
Afstotal = Afsfs + Agapfs;
Astr(1:Nfs,1:Nfs) = Afstotal;   % Add FS-FS connections; for now, count double connections as two edges
Astr(1:Nfs,Nfs+1:end) = Afsms;  % Add FS-MS connections

% get modularity!!!
[Sstr,U] = multileadevsplit(Astr);

% plot groups in 3D structure... and compute distance from centre
Dims = SIMPARAMS.net.PhysicalDimensions;
centre = Dims/2;
mediandist = zeros(thetaHngrps,1);
cdist = cell(thetaHngrps,1);
figure(106), clf,    
for loop = 1:thetaHngrps(bH)
    %inds = R(Sc==Igrp(loop));
    inds = RH(ScH==loop);
    xyz = [SIMPARAMS.net.MS.Position(inds,1),SIMPARAMS.net.MS.Position(inds,2),SIMPARAMS.net.MS.Position(inds,3)];
    cdiff = xyz - repmat(centre,thetaHgrpsizes{bH}(loop),1);
    cdist{loop} = sqrt(cdiff(:,1).^2 + cdiff(:,2).^2 + cdiff(:,3).^2);
    mediandist(loop) = median(cdist{loop});
    subplot(211),plot3(SIMPARAMS.net.MS.Position(inds,1),SIMPARAMS.net.MS.Position(inds,2),SIMPARAMS.net.MS.Position(inds,3),'.','Color',clrs(clrstep*loop,:))
    grid on, hold on
    subplot(212), hold on, [h,stats] = cdfplot(cdist{loop}); xlabel('x = distance from centre')
    set(h,'Color',clrs(clrstep*loop,:))
end  


save posthoc_analysis
