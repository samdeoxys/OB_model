ind = 1;
SI = masterdata(ind).InputCurrent;
SM = masterdata(ind).Mitral;
SGcd = masterdata(ind).GraDistal;
%%
j=357;
figure
subplot(3,1,1)
plot(SI.ImitgradistAMPA(j,:),'LineWidth',2)
hold on
plot(SI.ImitgradistNMDA(j,:),'LineWidth',2)
plot(SI.ImitgradistVDCC(j,:),'LineWidth',2)
%plot(SI.Igradistmit(j,:),'LineWidth',2)
hold off
legend('AMPA','NMDA','VDCC','GABA')
subplot(3,1,2)

i=357;
%plot(SM(i).V,'LineWidth',2)
hold on
plot(SGcd(i).V,'LineWidth',2)
hold off
%legend('Vmit','Vgcd')


subplot(3,1,3)
for i=357:362
%for i = 1:5 
    hold on
    plot(SGcd(i).V,'LineWidth',1.5)
end
hold off
legend
%%
populationmin = zeros(720,1);
for i=1:720
    populationmin(i) = min(SGcd(i).V(1,1000:end));
end
figure
histogram(populationmin)
title('V equilibrium for GCDs with changed Vrest')
mean(populationmin)
%%
temporalmean = zeros(720,1);
for i=1:720
    temporalmean(i) = mean(SGcd(i).V(1,1000:end));
end
figure
histogram(temporalmean)
%title('V temporal mean for GCDs with changed Vrest = -0.054')
title('V temporal mean for GCDs with decreased leaky = 0.62')
mean(temporalmean)
%% Vgcd pairwise correlation
Vgcd = get_neuron_field(SGcd,'V');
trim = 1000;
corrmat_Vgcd = corr(Vgcd(:,trim:end)');
corrmat_Vgcd = triu(corrmat_Vgcd);
corrvec_Vgcd = corrmat_Vgcd(:);
corrvec_Vgcd = corrvec_Vgcd(find(corrvec_Vgcd~=0&corrvec_Vgcd~=1));
nbins = 100;
figure
histogram(corrmat_Vgcd(:),nbins)
title('pairwise correlation of Vgcd, wmod = 0.62')
%% convolved Smit correlation and distance
taudist = 100; % unit = 0.1ms
[dist_mean,dist ]= get_dist(SM,taudist);
[corrmean,corrmat_Smit] = get_synchrony(SM,taudist,0);
%% mit firing properties
figure
Smit = get_neuron_field(SM,'S');
Smit_mean_vec = mean(Smit(:,trim:end),2)*10000;
Smit_std_vec = std(Smit(:,trim:end),0,2)*10000;
Smit_cv_vec = Smit_std_vec./Smit_mean_vec;
stat_Smit.mean = mean(Smit_mean_vec);
stat_Smit.std = mean(Smit_std_vec);
stat_Smit.cv = mean(Smit_cv_vec,'omitn');
Smit_master = {Smit_mean_vec,Smit_std_vec,Smit_cv_vec};
names = fieldnames(stat_Smit);
nbins = 15;
nfields = length(names);
for i=1:nfields
    subplot(nfields,1,i)
    histogram(Smit_master{i},nbins)
    title(['Smit ',names{i}])
    hold on
    line([stat_Smit.(names{i}),stat_Smit.(names{i})],get(gca,'YLim'),'color','r','LineS',':','LineW',1.0)
    legend('',num2str(stat_Smit.(names{i})))
    hold off
end
%%
ind = [1,5]
names = ["0.06","0.03"];
cmat = ones(3)*1/45; 
for i=1:2
    SI = masterdata(ind(i)).InputCurrent;
    SM = masterdata(ind(i)).Mitral;
    SGcd = masterdata(ind(i)).GraDistal;
    figure
    title(['Igaba, vrest = ',names(i)])
    hold on
    for j=1:45
        plot(SI.Igradistmit(j,:),'color',1-j*cmat(1,:));
    end
    hold off
end
%% animation of histogram 
SGcd = cell(1,2);
Vgcd = cell(1,2);
SMit = cell(1,2);
VMit = cell(1,2);
SI = cell(1,2);
plotax = cell(1,2);
numhighext=cell(1,2);


%figure
for i=1:2
    SGcd{i} = masterdata(i).GraDistal;
    SMit{i} = masterdata(i).Mitral;
    Vgcd{i} = get_neuron_field(SGcd{i},'V');
    Vmit{i} = get_neuron_field(SMit{i},'V');
    SI{i} = masterdata(i).InputCurrent;
    
    numhighext{i} = get_numhighext(SMit{i},SGcd{i});
    [~,indmax] = max(numhighext{i}(:,1));
    [~,indmin] = min(numhighext{i}(:,1));
    plot_mitral(SMit{i},SI{i},[indmax,indmin])
    
    
    %plotax{i} = subplot(2,1,i);
    %set(plotax{i},'XLim',[-20e-3,20e-3])
    %set(plotax{i},'XLim',[0,1])
    %hold on
    %xlim([-0.08,0.07])
    %xlim([-0.08,-0.01])
    %ylim([0,0.4])
end

%% high low compare prelease --- prototype for function high_low_compare
figure
hold on
plot(Pseperated{indmax,1})
plot(Pseperated{indmin,1})
hold off
plotbrowser on
datacursormode on
highgcind = find([SGcd{2}.Vrest]~=-0.074); 
lowgcind = find([SGcd{2}.Vrest]==-0.074);
mean_highgc_prelease = mean(SI{2}.Prelease(highgcind,:),1);
mean_lowgc_prelease =  mean(SI{2}.Prelease(lowgcind,:),1);
figure('un','norm','pos',[0.1,0.1,0.6,0.7])
hold on
title('mean gc groups prelease')
plot(mean_highgc_prelease,'LineW',1.5)
plot(mean_lowgc_prelease,'LineW',1.5)
hold off
legend('high','low')
plotbrowser on
datacursormode on

%%
for i=1:2
    [high,low]=high_low_compare(SGcd{i},SI{i},'CCaBase');
end


%%
for tt = 1:(tsim/dt)
    for i=2
        subplot(2,1,i)
        %plothist = histogram(Vgcd{i}(:,tt),'norm','prob','FaceC','b')
        plothist = histogram(Vmit{i}(:,tt),'norm','prob','FaceC','b');
        %plothist = histogram(SI{i}.Igradistmit(:,tt),'norm','prob','FaceColor','b');
        %plothist = histogram(SI{i}.Prelease(:,tt),'norm','prob','FaceColor','b');
        
        drawnow limitrate
        delete(plothist)
    end
end
   


        