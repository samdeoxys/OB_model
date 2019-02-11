function varargout = analyze_synchrony(mit,gc,taudist,fig_title)
% ---Sam Zheng---
%plots: histogram of equilibrium Vgc, with the mean; 
%       histogram of pairwise correlation of Vgcd, with the mean,  
%       histogram of pairwise correlation of convolved Smit, with the mean
%       histogram of pairwise distance of convolved Smit, with the mean
%input : 
% mit = struct array Mitral 
% gc = struct array GraDistal or GraProximal
% taudist = tau for calculating the distance function
% fig_title = title of the figure
%output: 
% out ---- a struct of the mean of the histograms
out = struct();
nmit = length(mit);
ngc = length(gc);
trim = 1000;
Vgc_populationmin = zeros(ngc,1);

figure('un','norm','Pos',[0.01,0.01,0.65,0.9])
subplot(2,2,1) %%% Vgc excitability
for i=1:ngc
    Vgc_populationmin(i) = min(gc(i).V(1,trim:end));
end
histogram(Vgc_populationmin,'norm','prob')
title(['V equilibrium for GCs ',fig_title])
Vgc_equilibrium_mean = mean(Vgc_populationmin);
out.Vgc_excitability = Vgc_equilibrium_mean;
out.Vgc_excitability_mat = Vgc_populationmin;
hold on
line([Vgc_equilibrium_mean,Vgc_equilibrium_mean],get(gca,'YLim'),'color','r','LineS',':','LineW',1.5)
legend('',num2str(Vgc_equilibrium_mean))
hold off
 

subplot(2,2,2) %%% Vgc correlation
Vgc = get_neuron_field(gc,'V');
corrmat_Vgc = corr(Vgc(:,trim:end)');
corrmat_Vgc = triu(corrmat_Vgc);
corrvec_Vgc = corrmat_Vgc(:);
corrvec_Vgc = corrvec_Vgc(find(corrvec_Vgc~=0&corrvec_Vgc~=1));
nbins = 100;
histogram(corrmat_Vgc(:),nbins,'norm','prob')
title(['pairwise correlation of Vgc, ',fig_title])
corr_Vgc_mean = mean(corrvec_Vgc);
hold on
line([corr_Vgc_mean,corr_Vgc_mean],get(gca,'YLim'),'color','r','LineS',':','LineW',1.5)
legend('',num2str(corr_Vgc_mean))
hold off
out.corr_Vgc = corr_Vgc_mean;


subplot(2,2,3) %%% mit distance
[dist_mean,dist]= get_dist(mit,taudist);
out.dist_Smit = dist_mean;
histogram(dist,nbins,'norm','prob')
title(['distance of convolved Smit',fig_title])
hold on
line([dist_mean,dist_mean],get(gca,'YLim'),'color','r','LineS',':','LineW',1.5)
legend('',num2str(dist_mean))
hold off

subplot(2,2,4) %%% mit correlation
[corrmean,corrmat_Smit] = get_synchrony(mit,taudist,0);
out.corr_Smit = corrmean;
histogram(corrmat_Smit,nbins,'norm','prob')
title(['pairwise correlction of convlved Smit ',fig_title])
hold on
line([corrmean,corrmean],get(gca,'YLim'),'color','r','LineS',':','LineW',1.5)
legend('',num2str(corrmean))
hold off

tightfig
varargout{1} = out;