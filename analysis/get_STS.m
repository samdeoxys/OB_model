function [sts,varargout] = get_STS(S,binsize,dt,tsim)
% calculate the spike train synchrony index of the network 
% from Brunel & Wang (2003)
% Input:
% S -- a struct array of neurons that can spike
% binsize -- binsize for histcount in ms
% dt, T -- time step of integration and total simulation time
fr_mat_original = get_fr(S,binsize,dt,tsim);
trim = 100; % cut what's before 100ms
fr_mat = fr_mat_original(:,floor(trim/binsize):end);
%varargout{1} = fr_mat_original;
fr_of_each_cell = mean(fr_mat,2);
%varargout{2} = fr_of_each_cell;
network_activity_across_time = sum(fr_mat);
mean_fr_across_cells = mean(fr_of_each_cell);
[r,lags] = xcorr(network_activity_across_time,0);
sts = r/mean_fr_across_cells^2;
sts = r/(mean(network_activity_across_time))^2;




