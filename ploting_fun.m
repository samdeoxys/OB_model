figure
hold on
plot(GraDistal(1).V)
%plot(GraProximal{1,1}.V)
plot(Mitral(1).V)
legend('gradistal.V','proximal.V','Mitral.V')
title('V')
hold off
figure
hold on
plot(InputCurrent.Igradistmit(1,:))
plot(InputCurrent.ImitgradistAMPA(1,:)+InputCurrent.ImitgradistNMDA(1,:)+InputCurrent.ImitgradistVDCC(1,:)-InputCurrent.IVDCCBase(1))

plot(InputCurrent.Ipcgraprox(1,:))
legend('ipsp on mitral', 'epsp on gcd','pcinput')
hold off
title('I')

figure
hold on 
plot(MitLFPs.VG);
plot(MitLFPs.extra);
legend('i lfp','e and i lfp')
title('LFP')
hold off

figure
hold on 
title('specific current')
plot(InputCurrent.ImitgradistAMPA(1,:))
plot(InputCurrent.ImitgradistNMDA(1,:))
plot(InputCurrent.ImitgradistVDCC(1,:)-InputCurrent.IVDCCBase(1))
plot(InputCurrent.ImitgraproxNMDA(1,:))
legend('ampa','nmda','vdcc','nmda-prox')

figure
hold on
title('activations')
plot(InputCurrent.Prelease(1,:));
plot(InputCurrent.CCa(1,:));
legend('Prelease','calcium')
hold off

firing_rate = sum(Mitral(1).S) / 0.7
%%
figure
hold on
for i=1:5
    plot(InputCurrent.ImitgraproxAMPA(i,:))
end
hold off

%% GraDs and grap connected to the mit
figure('Units','normalized','Position',[0.1,0.5,0.7,0.6])
hold on
plotbrowser('on')
colors = {'r','g','b'};
connected_graD = cell(1,3);
connected_graP = cell(1,3);
for i=1:3
    
    connected_graD{i} = find([masterdata(i).Mitral(1).ConnectionsGrad]);
    connected_graP{i} = find([masterdata(i).Mitral(1).ConnectionsGrap]);
    connected_graD_V = cell2mat(transpose({masterdata(i).GraDistal(connected_graD{i}).V}));
    connected_graP_V = cell2mat(transpose({masterdata(i).GraProximal(connected_graP{i}).V}));
    %connected_graP_V =
    %cell2mat(transpose({masterdata(i).GraProximal(connected_graP).V}));
    %%%uncomment this for future
    plot(mean(connected_graD_V),'LineStyle','--','color',colors{i},'LineWidth',1.5)
    plot(mean(connected_graP_V),'LineStyle','-.','color',colors{i},'LineWidth',1.5)
    plot(masterdata(i).InputCurrent.Igradistmit(1,:),'LineStyle','-','color',colors{i},'LineWidth',1.5)
    plot(masterdata(i).InputCurrent.Igraspmit(1,:),'LineStyle',':','color',colors{i},'LineWidth',1.5)
    
    mean_connected_graP_V = mean(connected_graP_V);
    
    ind = find(connected_graP_Stot{1,i}>0); %find the times when some gcp spiked
    maxnumsp = max(connected_graP_Stot{1,i}); %get the max num of gcp spiking together
    markersizes = 60+200*connected_graP_Stot{1,i}(ind)/maxnumsp; % set the marker sizes relative to its ratio to max number of gcp spiking together
    scatter(ind,mean_connected_graP_V(ind),markersizes,'filled');
    
end
legend
title('grad and grap connected to mit')
hold off

%% test for graP spike hypothesis
% the idea is to plot the time when some connected gcp spiked when plotting the EPSC; the marker
% size is determined by the number of cells fired at the time; but actually
% there's no spike; 
connected_graP_S = cell(1,3);
connected_graP_Stot = cell(1,3);
for i=1:3
    connected_graP_S{1,i} = cell2mat(transpose({masterdata(i).GraProximal(connected_graP{i}).S}));
    connected_graP_Stot{1,i} = sum(connected_graP_S{1,i});
end


%% current
figure('Units','normalized','Position',[0.1,0.5,0.7,0.6])
hold on
plotbrowser('on')
colors = {'r','g','b'};
markerstyle = {'o','+','d'};
for i=1:3
    ind = find(connected_graP_Stot{1,i}>0); %find the times when some gcp spiked
    maxnumsp = max(connected_graP_Stot{1,i}); %get the max num of gcp spiking together
    meanmitgradistAMPA = mean(masterdata(i).InputCurrent.ImitgradistAMPA);
    markersizes = 60+200*connected_graP_Stot{1,i}(ind)/maxnumsp; % set the marker sizes relative to its ratio to max number of gcp spiking together
    scatter(ind,meanmitgradistAMPA(ind),markersizes,'filled');
    
    plot(meanmitgradistAMPA,'LineWidth',2.0,'color',colors{i})
    plot(mean(masterdata(i).InputCurrent.Iext),'LineWidth',1.0,'LineStyle',':','color',colors{i})
    plot(masterdata(i).InputCurrent.ImitgradistAMPA(1,:),'LineWidth',1.0,'LineStyle',':','color',colors{i})
end
legend
title('MeanMitGradistAMPA + GraP spikes + Iext')
hold off
%% measurement the difference in each connection matrix: for the k mits connected to 1 gcd, sum the excitedness of the mits, plot the distribution
edges = 1:25;
for i=1:3
indicesum_connected_mit = cellfun(@sum,{masterdata(i).GraDistal.Connections});
figure
hold off
histogram(indicesum_connected_mit,edges,'Normalization','Probability')
ylim([0,0.16])
end

%% lfp
for i=1:3
    figure
    hold off
    plot(masterdata(i).MitLFPs.GradistMitGlobal)
end
%% population average firing rate
bin = 10;
figure 
hold on
for i=1:3
    mitrals_transposed = transpose({masterdata(i).Mitral.S})
    firing_mat = cell2mat(mitrals_transposed);
    FR(i,:) = sum(firing_mat,1);
    plot(smoothdata(FR(i,:)))
end
legend
title('population average firing rate')
%% mitral V
figure
hold on
colors = {'r','g','b'};

for i=1:3
    plot(masterdata(i).Mitral(1).V,'LineStyle','-','LineWidth',2,'color',colors{i})
    plot(masterdata(i).InputCurrent.Iext(1,:),'LineStyle',':','LineWidth',1.5,'color',colors{i})
end
legend
title('mitral 1')
hold off

%% grad and grap V
figure
hold on
colors = {'r','g','b'};

for i=1:3
    plot(masterdata(i).GraDistal(1).V,'LineStyle','-','LineWidth',2,'color',colors{i})
    plot(masterdata(i).GraProximal(1).V,'LineStyle',':','LineWidth',1.5,'color',colors{i})
end
legend
title('gra 1')
hold off

%% mit connected to connected gra (need the previous)
for i=1:3
    for j =1:length(connected_graD)
        ind = connected_graD(j);
        mit_connected_to_connected_graD(i,j) = sum(sum(find(masterdata(i).GraDistal(ind).Connections)));
    end
    for j=1:length(connected_graP)
        ind = connected_graP(j);
        mit_connected_to_connected_graP(i,j) = sum(sum(find(masterdata(i).GraProximal(ind).Connectionsmit))); %%remember to change to GraProximal
    end
end
mean(mit_connected_to_connected_graD,2)
mean(mit_connected_to_connected_graP,2)
%% test for synchrony hypothesis: 1 spike field coherence
params.tapers = [5 9];
params.fpass=[3 100];
params.Fs=10000;
params.err = [2 1]; %only work if err(1)=2, dont know why
params.trialave=0;
fscorr = 1;
C = cell(3,param.nMitral);
trim = 100;
figure
xlabel('Frequency')
ylabel('Coherence')
hold on
legendnames = {'mit.no1.data1','mit.no1.data2','mit.no1.data3'};
ntp = length(masterdata(1).MitLFPs.VG);
dataV = zeros(ntp-100-trim+1,1);
for i=1:3
    for j = 1:param.nMitral
        dataV = detrend(masterdata(i).MitLFPs.VG(:,trim:ntp-100))';
        [C{i,j},phi,S12,S1,S2,f,zerosp,confC,phistd,~]=coherencypb(dataV,masterdata(i).Mitral(j).S(trim:end-100)',params,fscorr);
        %[C{i,j},phi,S12,S1,S2,f,zerosp,confC,phistd,~]=coherencypt(masterdata(i).MitLFPs.VG,find(masterdata(i).Mitral(j).S)/10000,params,fscorr);
    end
    plot(f,C{i,1})
    
    legend(legendnames{1:i})
end
hold off
figure 
hold on
xlabel('Frequency')
ylabel('Coherence')
title('Average SFC')
for i=1:3
    plot(f,mean(cell2mat(C(i,:)),2))
    legend(legendnames{1:i})
end
hold off
%% what proportion of mit and grap fired together in a 10 dtpt window
figure('Unit','Normalized','Position',[0.1,0.5,0.7,0.6])
hold on 
datacursormode on
plotbrowser('on')
window = 14;
nbins = ntp/window;
syncMitsp = zeros(3,nbins-1);
syncGraPsp = zeros(3,nbins-1);
mitspikemat = cell(3,1);
grapspikemat = cell(3,1);
Legend = {'mitspike1','grapspike1','mitspike2','grapspike2','mitspike3','grapspike3'};
for i=1:3
    mitspikemat{i} = cell2mat({masterdata(i).Mitral.S}');
    grapspikemat{i} = cell2mat({masterdata(i).GraProximal.S}');
    for j = 1:nbins-1
        syncMitsp(i,j) = length(find(sum(mitspikemat{i}(:,j:j+window-1),2)))/45;
        syncGraPsp(i,j) = length(find(sum(grapspikemat{i}(:,j:j+window-1),2)))/720;
    end
    plot(syncMitsp(i,:),'LineWidth',1.5,'color',colors{i},'LineStyle','-')
    plot(syncGraPsp(i,:),'LineWidth',1.5,'color',colors{i},'LineStyle',':')
    legend(Legend{1:2*i})
end
title('proportion of mit/grap that fired together in a 10dtpt window')
hold off


%% plot several mitV
ind = [1,5,10,15,20,25,30,35,40];
figure
plotbrowser on
for i=1:length(ind)
    plot(masterdata(1).Mitral(ind(i)).V,'LineWidth',1)
    hold on
end

datacursormode on
legend
hold off
%% plotting prelease and calcium and VgraP 3
figure
hold on
plotbrowser on
for i = 1:3
     plot(masterdata(i).GraProximal(3).V)
     plot(masterdata(i).InputCurrent.Prelease(3,:))
     plot(masterdata(i).InputCurrent.CCa(3,:))
     plot(masterdata(i).InputCurrent.ImitgraproxNMDA(3,:))
     plot(masterdata(i).InputCurrent.ImitgraproxAMPA(3,:))
     plot(masterdata(i).InputCurrent.ImitgradistVDCC(3,:))
%     plot(masterdata(i).InputCurrent.mgradist(3,:))
%     plot(masterdata(i).InputCurrent.hgradist(3,:))
    plot(masterdata(i).InputCurrent.mgradist(3,:).*masterdata(i).InputCurrent.hgradist(3,:))
     plot(masterdata(i).Mitral(1).V)
     plot(masterdata(i).InputCurrent.Igradistmit(1,:))
end
legend
legend({'Vprox','P','Ca','mitgraproxNMDA','mitgraproxAMPA','mitgradistVDCC','m*h','Vmit','gaba'})
%% plot Ipcgraprox
figure
hold on
plotbrowser on
for i=1:3
    plot(masterdata(i).InputCurrent.Ipcgraprox(3,:),'color',colors{i})
end
legend

%% plot connected mit spike before the second spike of graprox3
figure
hold on
plotbrowser on
range = 500:1500;
for i = 1:3
    ind = find(masterdata(i).GraProximal(3).Connectionsmit);
    totsp = sum(mitspikemat{1}(ind,range));
    totsp = totsp>0;
        scatter(range,i*totsp,colors{i})
end
title('spikes of mit connected to graprox before its second spike')
%% average distance between two mits
distances = zeros(1,3);
allcomb = combnk(1:length(masterdata(1).Mitral),2);
%allcombmask = rand(length(allcomb),1);
%allcomb = allcomb(allcombmask <0.01,:);
numcomb = size(allcomb,1);
taudist = 30;
grapspikemat = cell(3,1);
for i=1:3
    grapspikemat{i} = cell2mat({masterdata(i).GraDistal.S}');
    distemp = 0;
    for j =1:numcomb
        row1 = allcomb(j,1);
        row2 = allcomb(j,2);
        distemp = distemp + spike_train_synchrony(mitspikemat{i}(row1,:),mitspikemat{i}(row2,:),taudist);
    end
    distances(i) = distemp / numcomb;
end
%% noise
figure
hold on
plotbrowser on
for i=1:3
    %plot(sum(masterdata(i).InputCurrent.Vnoisegraprox))
    plot(sum(masterdata(i).InputCurrent.Vnoisegradist))
end
legend
%%
figure
hold on
plot(Mitral(1).V)
plot(GraProximal(1).V)
legend
plotbrowser on
hold off

%% 
figure
hold on
for i=1:length(masterdata)
    plot(masterdata(i).Mitral(1).V)
    plot(masterdata(i).InputCurrent.ImitgradistAMPA(1,:))
    plot(masterdata(i).InputCurrent.ImitgradistNMDA(1,:))
    plot(masterdata(i).InputCurrent.ImitgradistVDCC(1,:))
    plot(masterdata(i).InputCurrent.Igradistmit(1,:))
    plot(masterdata(i).GraDistal(1).V)
    plot(masterdata(i).InputCurrent.Prelease(1,:))
end
hold off
legend('Vmit','AMPA','NMDA','VDCC','gaba','Vgradist','prelease','Vmit','AMPA','NMDA','VDCC','gaba','Vgradist','prelease','Vmit','AMPA','NMDA','VDCC','gaba','Vgradist','prelease')
%% 
