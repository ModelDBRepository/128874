function net = BuildStriatumNetwork(net,physiology,sim)

net.Cctms = int32(0:net.MS.N-1)';
net.Cctms_b = int32(0:net.MS.N)';
net.Cctms_d = int32(zeros(net.MS.N,1));
net.Cctms_w = ones(length(net.Cctms),1) .* 6.1;

net.Cctfs = int32(0:net.FS.N-1)';
net.Cctfs_b = int32(0:net.FS.N)';
net.Cctfs_d = int32(zeros(net.FS.N,1));
net.Cctfs_w = ones(length(net.Cctfs),1) .* 6.1;

% -------------------------------------------------------------------------
% connect the MS collaterals
probfile = 'MSNMSN_interconnection_data_v2.mat';
[net.Cmsms net.Cmsms_b net.Cmsms_d net.Cmsms_w] = physical_connections(net.MS.Position, net.MS.Position, probfile, net.MS.N, net.MS.N, physiology.MS_MS.maxdelay, physiology.MS_MS.baseweight, physiology.MS_MS.baseweightSD, net.PhysicalDimensions, sim.dt);

% connect the FS -> MS collaterals
probfile = 'FSMSN_interconnection_data_v2.mat';
[net.Cfsms net.Cfsms_b net.Cfsms_d net.Cfsms_w] = physical_connections(net.FS.Position, net.MS.Position, probfile, net.FS.N, net.MS.N, physiology.FS_MS.maxdelay, physiology.FS_MS.baseweight, physiology.FS_MS.baseweightSD, net.PhysicalDimensions, sim.dt);

% connect the FS collaterals
probfile = 'FSFS_interconnection_data_v2.mat';
[net.Cfsfs net.Cfsfs_b net.Cfsfs_d net.Cfsfs_w] = physical_connections(net.FS.Position, net.FS.Position, probfile, net.FS.N, net.FS.N, physiology.FS_FS.maxdelay, physiology.FS_FS.baseweight, physiology.FS_FS.baseweightSD, net.PhysicalDimensions, sim.dt);

% connect the FS gap junctions
probfile = 'FSgap_interconnection_data_v2.mat';
[net.Pgapfs net.Cgapfs net.Cgapfs_b net.Cgapfs_w] = physicalGapJunctions(net, physiology, sim, probfile);

% -------------------------------------------------------------------------
function [C C_b C_d C_w] = physical_connections(sourcecoords, targetcoords, probfile, N_source, N_target, maxdelay, baseweight, baseweightSD, dims, dt)

% load the probability function and parameters
load(probfile); 

% get the total number of terminals for each cell, and the total number of
% synapses
bounds = ones(N_source+1,1);
targetcells = [];
weight = [];
distance = [];

for i = 1:N_source
%     disp(N_source-i)
    % get the distance between the source cell and all the target cells
    d = get_distance(sourcecoords(i,:),targetcoords);

    % get the target cells and set some random weights
    [t w] = prob(d, i, N_target, exp_best_fun, exp_best_coeffs, baseweight ,baseweightSD, dims, sourcecoords(i,:));
    
    % get the bounds
    bounds(i+1) = bounds(i) + length(w);
    
    % fill the target cell and weight vectors
    targetcells = [targetcells; t];
    weight = [weight; w];

    % save the distance between the source and target cells. Needed to get
    % the delays as a function of distance.
    distance = [distance; d(targetcells(bounds(i):bounds(i+1)-1))];
end

% trim the vectors
% targetcells = targetcells(targetcells>0);
% weight = weight(weight>0);
% distance = distance(distance>0);

% get the transmission delay as a function of distance between the source
% and target cells
delay = round((distance ./ max(distance) .* maxdelay) / dt)+1;

% set the output variables
C = int32(targetcells-1);
C_b = int32(bounds-1);
C_d = int32(delay);
C_w = weight;

% -------------------------------------------------------------------------
function [t w] = prob(d, src, N_target, exp_best_fun, exp_best_coeffs, baseweight, baseweightSD, dims, sourcecoords)
% calculate the probability of making a connection, as a function of
% distance
p =  exp_best_fun(exp_best_coeffs,d);
p(src) = 0;     % no self-connection
t = find(rand(numel(p),1) < p);
w = (baseweight + randn(length(t),1) .* baseweightSD);

% make the connections
% count = 0;
% t_ = zeros(N_target,1);
% for i = 1:round(N_target) % go round extra time to account for edges
%     targcell = round((rand*N_target)+0.5); % a random number between 1 and N (number of cells)
%     if rand < p(targcell)
%         count = count + 1;
%         t_(targcell) = t_(targcell)+1;
%     end
% end
% x = sum(t_);
% 
% % calculate the weight, taking into account the number of synapses made
% T = [(1:N_target)' t_]; % gives the cell number and number of synapses
% t = T((T(:,2)>0),1);    % get the cells that have recived synapses
% w = (baseweight + randn(length(t),1) .* baseweightSD) .* T(t,2);

% -------------------------------------------------------------------------
function [Pgapfs Cgapfs Cgapfs_b Cgapfs_w] = physicalGapJunctions(net, physiology, sim, probfile)
load(probfile);

% get all possiable pairs of cells
[y d Nt] = posspermutations(net.FS.N, net.FS.Position);

% get the probabiity of a gap junction for each pair. Note, this is an
% inline function loaded from the probfile
p =  exp_best_fun(exp_best_coeffs,d); 
r = rand(length(p),1);
inds = find(r<(p.*1));
x = y(inds,:);

% find the number of gap juntions for each cell
numcon = zeros(net.FS.N,2);
for i = 1:net.FS.N
    numcon(i,:) = [i length(find(x==i))];
end

[Ngaps,jnk] = size(x);
Pgapfs = x;
x = [];

weights = physiology.FS_gap.baseweight + randn(Ngaps,1) .* physiology.FS_gap.baseweightSD; 

x = [Pgapfs(:,1) (1:Ngaps)' weights; Pgapfs(:,2) (1:Ngaps)' weights];

[jnk,i] = sort(x(:,1)); 
x = x(i,:);

Cgapfs = x(:,2);
Cgapfs_w = x(:,3);

% get the bounds
Cgapfs_b = zeros(net.FS.N,1);
for i = 1:net.FS.N
    gind = find(x(:,1) == i, 1, 'first');
    if isempty(gind)
        if i == 1
            Cgapfs_b(i) = 1;
        else
            Cgapfs_b(i) = Cgapfs_b(i-1);
        end
    else
        Cgapfs_b(i) = find(x(:,1) == i, 1, 'first');
    end
end

Pgapfs = int32(Pgapfs-1);
Cgapfs = int32(Cgapfs-1);
Cgapfs_b = int32([Cgapfs_b; length(x)+1])-1;

% -------------------------------------------------------------------------
function [y d Nt] = posspermutations(N, coords)
maxpercell = N-1;               % Maximum number of connections per cell
Nt = ((N.* maxpercell) ./2);    % total number of possiable gap junctions in the network

% y is the list of all possiable pairs
y = zeros(Nt,2);
count = 1;
for i = 1:N
    for j = i+1:N
        x = [i, j];
        y(count,:) = [min(x) max(x)];
        count = count + 1;
    end
end

[jnk,i] = sort(y(:,1));
y = y(i,:);

% get the distance between the pairs
d = zeros(Nt,1);
for i = 1:Nt
    a = coords(y(i,1),1) - coords(y(i,2),1);
    b = coords(y(i,1),2) - coords(y(i,2),2);
    c = coords(y(i,1),3) - coords(y(i,2),3);
    d(i) = sqrt(a.^2 + b.^2 + c.^2);
end

% -------------------------------------------------------------------------
% get the distance between the source and target cells
function d = get_distance(sourcecoords,targetcoords)
[n,m] = size(targetcoords);
a = targetcoords(:,1) - sourcecoords(1);
b = targetcoords(:,2) - sourcecoords(2);
c = targetcoords(:,3) - sourcecoords(3);
d = sqrt(a.^2 + b.^2 + c.^2);
