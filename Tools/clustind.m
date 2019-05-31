function [gamma,gammaG] = clustind(CIJ);

% input:  CIJ    = connection/adjacency matrix
% output: gamma  = cluster index for each vertex
%         gammaG = cluster index for entire graph

% Compute cluster index.
% Watts/Strogatz use un-directed graphs, here we use directed graphs.
% Find the immediate neighbors of a given vertex (both in and out), 
% then determine how many connections exist between them out 
% of all possible.  Note: immediate neighbors are all those vertices that
% either send or receive an edge from the target vertex.
% Written by Olaf Sporns, Indiana University, 2002
% Bug fixed (handling orphan neigbors, previously gamma = NaN) - 2004
% - Thanks Constantin

N = size(CIJ,1);

gamma = [];
% loop over all vertices
for v=1:N
   [nb] = find(CIJ(v,:) + CIJ(:,v)');
   if (~isempty(nb))
      gamma = [gamma sum(sum(CIJ(nb,nb)))./(length(nb)^2-length(nb))];
   end;
end;

% handle nodes that are connected to only a single other node
indices = find(isnan(gamma));
gamma(indices) = 0;

gammaG = mean(gamma);
