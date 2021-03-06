function fr_mat = get_fr(S,varargin)
%%% ----Sam Zheng---
%%% calculate the PSTH matrix for a struct array of neurons that can spike
%%% input:
%%% S -- a struct array of neurons that can spike
%%% (optional) edges -- edges for histcounts
%%% output: 
%%% fr_mat -- psth matrix, ncells * nbins
if exist(varargin)
    edges = varargin{1}
else
    edges = 1:100:7001; %data pts, =10ms
end

Smat = get_neuron_field(S,'S');
fr_mat = zeros(length(S),length(edges)-1);
for i = 1:length(S)
    inds = find(S(i).S);
    fr_mat(i,:) = histcounts(inds,edges)
    fr_mat(i,:) = fr_mat(i,:) / (edges(i+1)-edges(i)) * 10000;
end
