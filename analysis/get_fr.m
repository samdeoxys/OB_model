function fr_mat = get_fr(S,varargin)
%%% ----Sam Zheng---
%%% calculate the PSTH matrix for a struct array of neurons that can spike
%%% input:
%%% S -- a struct array of neurons that can spike
%%% (optional) binsize -- binsize for histcount in ms
%%% (optional) dt, T -- time step of integration and total simulation time
%%% output: 
%%% fr_mat -- psth matrix, ncells * nbins
if nargin > 1
    binsize = varargin{1};
    dt =varargin{2};
    tsim = varargin{3};
    edges = 1:floor(binsize/dt):floor(tsim/dt);
    
else
    edges = 1:100:7001; %data pts, =10ms
    binsize = 10;
end

Smat = get_neuron_field(S,'S');
fr_mat = zeros(length(S),length(edges)-1);
for i = 1:length(S)
    inds = find(S(i).S);
    fr_mat(i,:) = histcounts(inds,edges);
    fr_mat(i,:) = fr_mat(i,:) * (1000/binsize);
end
