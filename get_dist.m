function [dist_mean,varargout]=get_dist(N,taudist)
%%% ---Sam Zheng---
%%% Convolve spike trains of a struct array of neurons with an alpha
%%% function and compute the pairwise distances, and the mean of the
%%% distances.
%%% Input: 
%%% N = struct array of neurons that can spike
%%% taudist = parameter for the convolution kernel
%%% Output:
%%% dist_mean = the mean of the pairwise distances
%%% (optional) dist = the vector of all the pairwise distances
allcomb = combnk(1:length(N),2);
numcomb = size(allcomb,1);
dist = zeros(1,numcomb);
for j =1:numcomb
    row1 = allcomb(j,1);
    row2 = allcomb(j,2);
    dist(j) = spike_train_synchrony(N(row1).S(1,:),N(row2).S(1,:),taudist);
end
dist_mean = mean(dist);
varargout{1} = dist;
end

function distance = spike_train_synchrony(spike1,spike2,tau) 
% calculate the van Rossum distance between two spike trains (0/1)
clear t
datarange = 1:10; 
K = vRkernel(datarange,tau);
x = conv(K,spike1);
y = conv(K,spike2);
distance = 1/tau * sum((x-y).*(x-y));
end

function K = vRkernel(t,tau)
K = heaviside(t).*exp(-t./tau);

end