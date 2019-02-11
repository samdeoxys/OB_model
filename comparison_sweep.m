%%% --Sam Zheng--
%%% For comparing the sweep results for two different states, specified by
%%% the struct switchslider. eg. PCinputOn = true vs false

clear Mitral GraProximal GraDistal param MitLFPs GraDistLFPs masterdata
input_file = 'OB_params_GCE_Fig4.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

%%


tauvdccslider = 18;%*linspace(0.25,2,9);%[18];
nmdaslider = 0.045;%0.04*linspace(0.5,3,9); % 0.02
gcdtauslider = 4;%linspace(1,10,5);
mittauslider = 7;%linspace(2,12,6);
% nmdatauslider = [4]; % 2 has wide pattern if wnmda=0.04
nmdatau2slider =75;%75*linspace(2/3,2,9); %75
ccathslider = linspace(1,2,9);


%sine
%pcparam1slider = 0.3;%linspace(10,100,5); %amplitude in the case of sine
%pcparam2slider = [0.015];%30;%15;%linspace(0.005,0.05,10); %synaptic weight, controlling baseline strength
%pcparam3slider = linspace(45,65,3);

%poisson
pcparam1slider = 90; 
pcparam2slider = [15]; %since 15 mit onto 1 gra
pcparam3slider = [0.8];
pcparam4slider = [5];

pctypeslider = ["poisson"];
etaslider = 0.28;%linspace(0.04,0.28,9);% linspace(0.05,0.25,5);
Vgcdslider = -0.074;
Vgcpslider = -0.074;


switchslider.value = [true,false];
switchslider.name = 'PCinputON';
slider.value = ccathslider;
slider.name = 'ccath';
    
T = table();
subT = table();
masterdata(1:length(switchslider.value),1:length(slider.value)) = struct();
sweepfig = figure('un','N','pos',[0.1,0.1,0.4,0.7]);
subplot1(2,1,'Gap',[.01,.01],'XTickL','Margin','YTickL','Margin')


for i=1:length(switchslider.value)
    
    %create figure for spectrum, subplots and tighten them
    spectrumfig = figure('un','N','pos',[0.01,0.01,0.9,0.9]);
    subplot1(3,3,'Gap',[.01,.01])
    
    for j=1:length(slider.value)
        [Mitral GraProximal GraDistal param] = OB_network_GCE(input_file);
        param.(switchslider.name) = switchslider.value(i);
        
        [GraDistal.tauVDCC] = deal(tauvdccslider(1));
        %              [GraProximal.tau] =deal(gcptauslider(1));
             [Mitral.tau] = deal(mittauslider(1));
             [GraDistal.tau] = deal(gcdtauslider(1));
              [GraDistal.wNMDAMI] = deal(nmdaslider(1));
% %             % [GraDistal.wGABAGRSP] = deal(0);
%               [GraDistal.wAMPAMI] = deal(ampaslider(1));
%               [GraProximal.wAMPAMI] = deal(ampaslider(1));
               [GraProximal.wNMDAMI] = deal(nmdaslider(1));
%              [GraDistal.wVDCCMI] = deal(vdccslider(1));
%              [GraDistal.tauVDCC] = deal(tauvdccslider(1));
%               [Mitral.wGABAGR] = deal(gabaslider(1));
               [GraDistal.Vrest] = deal(Vgcdslider(1));
              [GraProximal.Vrest] = deal(Vgcpslider(1));
% %              %[GraProximal.FThresh] = deal(-0.05);
%               [GraDistal.tauAMPA1] = deal(ampatauslider(1));
%               [GraProximal.tauPROX1] = deal(1/2*ampatauslider(1));
%               [GraDistal.tauAMPA2] = deal(ampatau2slider(1));
%               [GraProximal.tauPROX2] = deal(2/3*ampatau2slider(1));
%               [GraDistal.tauNMDA1] = deal(nmdatauslider(1));
%               [GraProximal.tauNMDA1] = deal(nmdatauslider(1));
               [GraDistal.tauNMDA2] = deal(nmdatau2slider(1));
               [GraProximal.tauNMDA2] = deal(nmdatau2slider(1));
%               param.CChanceMitGraProx = deal(collateralslider(1));
              [GraDistal.CCaTh] = deal(ccathslider(j)); %1.5
%              [GraDistal.CCaR] = deal(0);
%              [GraDistal.tauCa] = deal(taucaslider(1));
% %             
%              param.rhoCaVDCC = deal(rhoCaNMDAslider(1));
%              param.rhoCaNMDA = deal(rhoCaVDCCslider(1));
% %             param.extparam = deal(extslider(1));
%              param.Mg_conc = 1;
              param.eta = etaslider(1);
%              param.gamma = 0.016;
%             param.NnoiseMit = 0.01*noiseslider(1);
%             param.NoiseGradist = 0.01*noiseslider(1);
%             param.NoiseGraProx = 0.01*noiseslider(1);
%             param.Inoise = 0.1*noiseslider(1);
            
 %           param.PCinputON = true;
            param.PCtype = pctypeslider(1);
            param.PCparam1 = pcparam1slider(1);
            param.PCparam2 = pcparam2slider(1);
            param.PCparam3 = pcparam3slider(1);
            param.PCparam4 = pcparam4slider(1);

        
        
        
        
        
        
        
        [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = ...
                        IandVLFP_GCE(Mitral,GraProximal,GraDistal,param);
        cmat = ones(3)*1/ngradist; % bw
        fs = 16;% fontsize

        trim = 1000; % trim beginning and end to avoid edge effects
            % trim:end-50
        L = length(timevec(trim:end-100));  % Length of simulation
        NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
        f = sampf/2*linspace(0,1,NFFT/2+1);

            %faxis = [50 90;20 60;10 40]; % frequency axis for power plots
        faxis = [10 90;10 90;10 90]; %now that paramslider is used, the axis has to encompass all

        mitFFTGI = fft(detrend(MitLFPs.GradistMitGlobal(trim:end-100),'constant'),NFFT)/L;
        mitFFT1I = fft(detrend(MitLFPs.GradistMit1(trim:end-100),'constant'),NFFT)/L;
        mitFFT2I = fft(detrend(MitLFPs.GradistMit2(trim:end-100),'constant'),NFFT)/L;
        mitFFT3I = fft(detrend(MitLFPs.GradistMit3(trim:end-100),'constant'),NFFT)/L;

        mitFFTEXTRA = fft(detrend(MitLFPs.extra(trim:end-100),'constant'),NFFT)/L; %[[Sam]]
        gradFFTGI = fft(detrend(GraDistLFPs.MitGradistGlobal(trim:end-100),'constant'),NFFT)/L; %[[Sam]]

        mitFFTGV = fft(detrend(MitLFPs.VG(trim:end-100),'constant'),NFFT)/L;
        mitFFT1V = fft(detrend(MitLFPs.V1(trim:end-100),'constant'),NFFT)/L;
        mitFFT2V = fft(detrend(MitLFPs.V2(trim:end-100),'constant'),NFFT)/L;
        mitFFT3V = fft(detrend(MitLFPs.V3(trim:end-100),'constant'),NFFT)/L;

        spectrumfig;
        subplot1(j)
        hold on        
        plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color','k','linestyle','--'); % bw
        plot(f,2*abs(gradFFTGI(1:NFFT/2+1)),'color','k','linestyle',':'); %[[Sam]]
        plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color','k','linestyle','-');
        %plot(f,1e-4*ones(1,length(f)),'k--')

        hold off
        xlim(faxis(1,:));ylim([0 2e-3])
        set(gca,'fontsize',fs)
                    %legend('ILFP - Global','VLFP - Global','ILFP - highest 15 MCs','ILFP - lowest 15 MCs')
        
        %xlabel('Frequency (Hz) ')
                    %title(['PCparam1 ',num2str(param.PCparam1),' Extparam',num2str(param.extparam)])
                   
        title(num2str(slider.value(j)))
        tableindex = (i-1)*length(slider.value) + j;
        T.(slider.name)(tableindex) = slider.value(j);
        [maxpower,ind] = max(2*abs(gradFFTGI(1:NFFT/2+1))); 
        T.f(tableindex) = f(ind);
        T.power(tableindex) = maxpower;
        T.(switchslider.name)(tableindex) = switchslider.value(i);
        
        
        %%
        masterdata(i,j).Mitral = Mitral;
        masterdata(i,j).GraDistal = GraDistal;
        masterdata(i,j).GraProximal = GraProximal;
        masterdata(i,j).InputCurrent = InputCurrent;
        masterdata(i,j).MitLFPs  = MitLFPs;
        masterdata(i,j).GraDistLFPs = GraDistLFPs;
         
                
    end
    spectrumfig %setting legend inside the loop causes trouble
    legend('ILFP','ELFP','VLFP') %[Sam]
    legend boxoff
    tightfig(spectrumfig)
    saveas(spectrumfig,[slider.name,'_',switchslider.name,num2str(switchslider.value(i)),'.jpg'])
    
    subT = T(T.(switchslider.name) == switchslider.value(i),:);
    sweepfig 
    subplot1(1)
    fignamef = [slider.name,' vs f.jpg'];
    hold on
    
    title(sweepfig.Children(1),fignamef) %interestingly if the plot is before the title it won't plot?
    plot(sweepfig.Children(1),subT.(slider.name),subT.f);
    hold off
    
    subplot1(2)
    hold on
    fignamep = [slider.name,' vs power.jpg'];
    title(sweepfig.Children(2),fignamep)
    plot(sweepfig.Children(2),subT.(slider.name),subT.power);
    hold off
end
leg = [switchslider.value(1),switchslider.value(2)];
% for i=1:2
%     legend(sweepfig.Children(i),string(switchslider.value))
% end



%% turning on datacursor for all figs
figHandles = findobj('type','figure');
for i=1:length(figHandles)
    datacursormode(figHandles(i))
end
fname = [slider.name,'_',switchslider.name];
savefig(figHandles,[fname,'_sweep'])
saveas(sweepfig,[fname,'_sweep.jpg'])
saveas(spectrumfig,[fname,'_spectrum.jpg'])