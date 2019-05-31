%%% script to analyse data-set using graph-cut method
clear all
load('..\RandomInput.mat'); load('..\RandomInput_SIMPARAMS.mat')

% binned spike count pars
binsize = [20:20:100, 200:100:1000];   %size of window in ms
%binsize = [10:10:100];

% correlation threshold
thetac = 0.3;   % for spike-count (words) -  has to be very low...

% Hamming threshold...
thetaH = 0.2;   % was 0.2

% limit on number of nodes required in graph for analysis
nlimit = 5;


% get indices of MS neurons
MSspks = out.STms; MSspks(:,1) = MSspks(:,1)+1; % change from zero-base to 1-base index 
idxs = unique(MSspks(:,1));
Nms = SIMPARAMS.net.MS.N;
Nmsidxs = numel(idxs);

% get indices of FS neurons
FSspks = out.STfs; FSspks(:,1) = FSspks(:,1)+1;
FSidxs = unique(FSspks(:,1));
Nfsidxs = numel(FSidxs);
Nfs = SIMPARAMS.net.FS.N;

% figure(1); clf; plot(MSspks(:,2), MSspks(:,1), '.')
% figure(10); clf; plot(FSspks(:,2), FSspks(:,1), '.')

% sim pars
simT = SIMPARAMS.sim.tfinal;
MSdt = SIMPARAMS.sim.dt;
simsteps = simT / MSdt;
simT_SI = simT * 1e-3;

% set up memory requirements
clear out                       % make some space!!

% binned rates
scount = cell(numel(binsize),1);
Sc_all = cell(numel(binsize),1);    % has to be cell array cos
Rn = cell(numel(binsize),1);        % retained neurons...
Ic_all = zeros(Nms,numel(binsize));
rmean = zeros(numel(binsize),1);        
rmax = zeros(numel(binsize),1); 
hmedian = zeros(numel(binsize),1);        
hmin = zeros(numel(binsize),1); 
hdiff = zeros(numel(binsize),1); 
Hamming_all = zeros(numel(binsize),1); 
Cmplx = zeros(numel(binsize),1);
grpsizes = cell(numel(binsize),1);
kgrpsizes = cell(numel(binsize),1);
dgr = cell(numel(binsize),1);
dgrRn = cell(numel(binsize),1);
Hamming = cell(numel(binsize),1);
kHamming = cell(numel(binsize),1);

for loop = 1:numel(binsize)
    loop
    bins = 0:binsize(loop):simT;
    scount{loop} = zeros(numel(bins),Nms);
    for j = 1:Nmsidxs
        currix = find(MSspks(:,1) == idxs(j));
        nspikes(j) = numel(currix);
        scount{loop}(:,idxs(j)) = histc(MSspks(currix,2),bins);
    end
    % analyse...
    %%%%%%%%%%%%%%%%%% analyse spike count (words) %%%%%%%%%%
    % do correlation on binary vectorised spike count
    scount{loop}(scount{loop} > 1) = 1;
    % count number of > 1...
    
    % Hamming distance (error)....
    HE = pdist(scount{loop}','hamming');
    hc = squareform(HE);    % turn it into matrix - with zeroes on diagonal
    Ac = (hc <= thetaH);    % threshold to get matrix..
    hmedian(loop) = median(HE(HE>0));    % store median of Hamming error that are not perfect
    hmin(loop) = min(HE(HE>0));  % store min Hamming error that are not perfect
    hdiff(loop) = hmedian(loop) - hmin(loop);    
    
    % correlation coefficient
%     rc = corrcoef(scount{loop});
%     rc(eye(size(rc))==1) = 0;  % no self-connections
%     rc(isnan(rc)) = 0;
%     Ac = (rc >= thetac);
%     rmean(loop) = sum(sum(rc.^2))./(Nms^2-Nms);    % store mean of r^2 
%     rmax(loop) = max(max(rc));  % store max r
   
    % strip out low linked nodes...
    dgr{loop} = sum(Ac);
    delN = find(dgr{loop} < 2);
    Rn{loop} = setdiff(idxs,delN);    %% index back into Rn to recover original index of MSN...
    Acs = Ac(Rn{loop},Rn{loop});
    dgrRn{loop} = sum(Acs);
    
    % do modularity on that graph
    edges(loop) = sum(sum(Acs)); nodes(loop) = numel(Rn{loop});
    if nodes(loop) > nlimit & edges(loop) > log(nodes(loop))   % giant component phase transition 
        % only group if there are enough nodes [how many???] and  edges to make it worthwhile...
        [Sc_all{loop},Uc] = multileadevsplit(Acs) ;% ,'refine'); % do refine step later when got time...
        ngrps(loop) = max(Sc_all{loop});
        siz = [];
        for i = 1:ngrps(loop)
            siz = [siz numel(find(Sc_all{loop} == i))];   
        end
        grpsizes{loop} = siz;
    else
        Sc_all{loop} = [];  % nothing to group!
        ngrps(loop) = 0;
        grpsizes{loop} = [];
    end
    
       
    % do k-means for comparison, using same number of groups as detected by
    % this method
    % use binary vectors with k-means as well, and compute Hamming....
    % [Ic_all(:,loop),Cc] = kmeans(scount{loop}',k,'Distance','Hamming'); 
    if ngrps(loop) == 0 
        k(loop) = 2; % graph-cut analysis found zero correlated groups, try k-means with just two!
    elseif numel(Rn{loop}) < Nmsidxs   
        k(loop) = ngrps(loop)+1; % do with one extra group for all discarded neurons?
    else
        k(loop) = ngrps(loop); % add neurons in graph-cut analysis
    end
    
    try 
        [Ic_all(:,loop),Cc] = kmeans(scount{loop}',k(loop),'Distance','Hamming');   % same number of groups...
        siz = [];
        for i = 1:k(loop)
            siz = [siz numel(find(Ic_all(:,loop) == i))];
            % and compute group Hamming
            cts = [scount{loop}(:,Ic_all(:,loop)==i)];
            kHamming{loop}(i) = median(pdist(cts','hamming'));
        end
        kgrpsizes{loop} = siz;
    catch
        % will throw error if finds empty cluster on first iteration  
        k(loop) = 0;
        kHamming{loop} = 0;
        kgrpsizes{loop} = 0;
    end

    
    %%% Hamming distance within groups vs whole retained network for
    %%% modularity method
    Hamming_all(loop) = median(pdist(scount{loop}(:,Rn{loop})','hamming'));

    for i=1:ngrps(loop)
        cts = [scount{loop}(:,Rn{loop}(Sc_all{loop}==i))];
        Hamming{loop}(i) = median(pdist(cts','hamming'));
    end
    
    if ngrps(loop)
        % compute "complexity" measure: eq 19 in Humphries et al (2009)
        Cmplx(loop) = sum(grpsizes{loop})/Nms * ngrps(loop)  * hdiff(loop);
    end
end

%%---------------------------------------------------------- analyse output


%%% only count grps bigger than 1%??
grpthresh = round(Nms * 0.01);
ngrpsT = zeros(numel(ngrps),1);
for i = 1:numel(ngrps)
    ngrpsT(i) = sum(grpsizes{i} > grpthresh);
end

%%% trend in grps, mean r^2, and max r with binsize
figure(2), clf
subplot(411),plot(binsize,ngrps,'+-'), hold on
plot(binsize,ngrpsT,'k+-')
ylabel('# groups')
subplot(412), plot(binsize,Cmplx,'+-'),
ylabel('Complexity of structure')
if any(rmean)
    subplot(413),plot(binsize,rmean,'+-'),ylabel('Mean r^2')
    subplot(414),plot(binsize,rmax,'+-'),ylabel('Max r')
else
    subplot(413),plot(binsize,hdiff,'+-'),ylabel('H error difference')
    subplot(414),plot(binsize,hmin,'+-'),ylabel('Min Hamming error')    
end
xlabel('bin size (ms)')

%%%
%btest = find(ngrps == max(ngrps)); 
% btest = find(ngrpsT == max(ngrpsT)); 
btest = find(Cmplx == max(Cmplx)); % pick binsize by most complex

% in case there's more than one, pick with greatest difference in H error from set
% of binsizes that resulted in groups!
if numel(btest) > 1 
    ix = find(hdiff(btest) == max(hdiff(btest)));
    btest = btest(ix);
end

%btest = 8;
Sc = Sc_all{btest};
Ct = scount{btest};
R = Rn{btest};
Ic= Ic_all(:,btest);
Hgrps = Hamming{btest}';    % intra-group Hamming distance
kHgrps = kHamming{btest}';

%%% sort groups into Hamming error order 
[Hs,Igrp] = sort(Hgrps);
[kHs,kIgrp] = sort(kHgrps);
grp_sz_srt = grpsizes{btest}(Igrp);

%%% get group membership in cell array
grp_members = {};
nD1 = zeros(ngrps(btest),1); nD2 = zeros(ngrps(btest),1);
for loop = 1:ngrps(btest)
   grp_members{loop} = R(Sc==loop); 
   for i = 1:grpsizes{btest}(loop)
       if any(grp_members{loop}(i) == SIMPARAMS.net.MS.D1inds) 
           nD1(loop) = nD1(loop) + 1;
       elseif any(grp_members{loop}(i) == SIMPARAMS.net.MS.D2inds) 
           nD2(loop) = nD2(loop) + 1;
       else
           warning('Not D1 or D2')
       end
   end
end


%%%%%% Uncomment below to get plots in Hamming order %%%%%%%%%%%%%%%%%%%%%
%%% plot out bins
bins = 0:binsize(btest):simT;
scts = zeros(numel(bins),ngrps(btest));
figure(3), clf
for loop = 1:ngrps(btest)
    %cts = [Ct(:,R(Sc==Igrp(loop)))]';
    cts = [Ct(:,R(Sc==loop))]';
    scts(:,loop) = [sum(cts)./grpsizes{btest}(loop)]';
    subplot(max(Sc),1,loop), bar(bins,sum(cts)./grpsizes{btest}(loop),'histc')
    title('Summed spike count for graph-method')
end

sctsK = zeros(numel(bins),k(btest));
figure(4), clf
for loop = 1:k(btest)
    % cts = [Ct(:,Ic == kIgrp(loop))]';
    cts = [Ct(:,Ic == loop)]';
    sctsK(:,loop) = [sum(cts)]';
    subplot(k(btest),1,loop), bar(bins,sum(cts),'histc')
    title('Summed spike count for k-means')
end

figure(5), clf, hold on, title('Grouped raster by graph-method')
% full colour version
clrs = colormap; [rw cl] = size(clrs); clrstep = floor(rw / ngrps(btest));
% alternating black and grey version
% clrs = repmat([0 0 0; 0.5 0.5 0.5],ngrps(btest),1); clrstep = 1;

ctr = numel(R);
for loop = 1:ngrps(btest)
    %inds = R(Sc==Igrp(loop));
    inds = R(Sc==loop);
    stgrp = []; ixgrp=[];
    for i = 1:numel(inds)
        indst = find(MSspks(:,1) == inds(i));
        st = MSspks(indst,2);
        stgrp = [stgrp; st];
        ixgrp = [ixgrp; ones(length(st),1)+ctr];
        ctr = ctr-1;
    end
    plot(stgrp, ixgrp,'.','Color',clrs(clrstep*loop,:));
    %fname1 = ['stgrp' num2str(Igrp(loop)) '.txt']; fname2 = ['ixgrp' num2str(Igrp(loop)) '.txt']; 
    %save(fname1,'stgrp','-ascii'); save(fname2,'ixgrp','-ascii')
end
% axis off


% plot groups in 3D structure...
figure(6), clf,    
for loop = 1:ngrps(btest)
    %inds = R(Sc==Igrp(loop));
    inds = R(Sc==loop);
    plot3(SIMPARAMS.net.MS.Position(inds,1),SIMPARAMS.net.MS.Position(inds,2),SIMPARAMS.net.MS.Position(inds,3),'.','Color',clrs(clrstep*loop,:))
    hold on
end  

% do raster based on k-means...
figure(7), clf, hold on, title('Grouped raster by k-means (group 1 at the top)')
clrs = hsv; [rw cl] = size(clrs); clrstep = floor(rw / k(btest));
ctr = numel(Ic);
for loop = 1:k(btest)
    %inds = find(Ic==kIgrp(loop));
    inds = find(Ic==loop);
    st = [];
    for i = 1:numel(inds)
        indst = find(MSspks(:,1) == inds(i));
        st = MSspks(indst,:);
        plot(st(:,2), ones(length(st),1)+ctr,'.','Color',clrs(clrstep*loop,:));
        % plot(st(:,2), ones(length(st),1)+ctr,'.');
        ctr = ctr-1;
    end
end

%%% save data
% fname = ['theta' num2str(thetac) '.mat'];
% save(fname)
save first_stage_analysis