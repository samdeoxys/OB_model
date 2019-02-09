clear Mitral GraProximal GraDistal param MitLFPs GraDistLFPs
input_file = 'OB_params_GCE_Fig4.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);
%[Mitral GraProximal GraDistal param] = OB_network_GCE(input_file);


%% construct sliders
ampaslider = 0.024; %linspace(0.01,0.082,4);%linspace(0.01,0.082,10); %0.024
nmdaslider = 0.02;%linspace(0.02,0.08,10);%linspace(0.02,0.08,10); %[0.04]; % 0.02
nmdatauslider = 4;%linspace(2,20,4);%[4]; % 2 has wide pattern if wnmda=0.04
nmdatau2slider =75; %75
ampatauslider = 2;%*linspace(0.1,1.4,20);%2
ampatau2slider = 3;%*linspace(0.7,3,20);  %3
%tau2slider = 2*tauslider;
gabaslider = linspace(0.002,0.02,10); % 0.008
mittauslider = 5;%linspace(1,10,4); %linspace(1,10,20); linspace(1,10,10);
gcdtauslider = 5;%linspace(1,10,10);%linspace(1,10,10);%5; %linspace(1,10,20);
gcptauslider = 5;%linspace(1,10,20);
ccathslider = 1.5; %1.5
extslider = [1.25];
vdccslider = 200;%linspace(0,400,20); %200
tauvdccslider = 18;%18*linspace(0.1,0.5,10);%18*linspace(0.25,2,20);%[18];
collateralslider = [0.2];
noiseslider = [0];
taucaslider = 5*linspace(0.25,2.5,4);%5;%5*linspace(0.25,2.5,20);%5;
rhoCaNMDAslider = 50;%50 * linspace(1/2,4,20);
rhoCaVDCCslider = 100;%50 * linspace(1/2,4,20); % 100



% [GraDistal.wAMPAMI] = deal(0.03);
% %[GraDistal.wNMDAMI] = deal(0.04);
% [GraProximal.wAMPAMI] = deal(0.04);
% [GraProximal.wNMDAMI] = deal(0.03);
slider2 = [];
pcparam2slider = [0];
pctypeslider = ["constant"];
pcparam1slider = [1];
% slider that needs to be changed whenever looping on a different variable
slider.value = cell(2,1);
slider.name = cell(2,1);
slider.value{2} = gabaslider;
slider.name{2} = 'wgaba';
slider.value{1} = taucaslider;
slider.name{1} = 'tauca';
%% run
masterdata = struct;
T = struct();
fieldnames = cell(length(slider.value{1}),1);
for i=1:length(slider.value{1}) % slider1 is the high mid low condition, slider2 is the sweep condition
    fieldnames{i} = [slider.name{1},'_',num2str(slider.value{1}(i))];
    fieldnames{i} = replace(fieldnames{i},'.','dot'); % fieldnames cannot contain dot!!!
    T.(fieldnames{i}) = table();
    for j=1:length(slider.value{2})
    
  %  if ~isempty(slider2)
   %     for j=1:length(slider2)
            [Mitral GraProximal GraDistal param] = OB_network_GCE(input_file);
             [GraProximal.tau] =deal(gcptauslider(1));
            [Mitral.tau] = deal(mittauslider(1));
            [GraDistal.tau] = deal(gcdtauslider(1));
              [GraDistal.wNMDAMI] = deal(nmdaslider(1));
%             % [GraDistal.wGABAGRSP] = deal(0);
              [GraDistal.wAMPAMI] = deal(ampaslider(1));
              [GraProximal.wAMPAMI] = deal(ampaslider(1));
              [GraProximal.wNMDAMI] = deal(nmdaslider(1));
             [GraDistal.wVDCCMI] = deal(vdccslider(1));
             [GraDistal.tauVDCC] = deal(tauvdccslider(1));
              [Mitral.wGABAGR] = deal(gabaslider(j));
              [GraDistal.Vrest] = deal(-0.074);
             [GraProximal.Vrest] = deal(-0.074);
%              %[GraProximal.FThresh] = deal(-0.05);
              [GraDistal.tauAMPA1] = deal(ampatauslider(1));
              [GraProximal.tauPROX1] = deal(1/2*ampatauslider(1));
              [GraDistal.tauAMPA2] = deal(ampatau2slider(1));
              [GraProximal.tauPROX2] = deal(2/3*ampatau2slider(1));
              [GraDistal.tauNMDA1] = deal(nmdatauslider(1));
              [GraProximal.tauNMDA1] = deal(nmdatauslider(1));
              [GraDistal.tauNMDA2] = deal(nmdatau2slider(1));
              [GraProximal.tauNMDA2] = deal(nmdatau2slider(1));
              param.CChanceMitGraProx = deal(collateralslider(1));
             [GraDistal.CCaTh] = deal(ccathslider(1)); %1.5
             [GraDistal.CCaR] = deal(0);
             [GraDistal.tauCa] = deal(taucaslider(i));
%             
             param.rhoCaVDCC = deal(rhoCaNMDAslider(1));
             param.rhoCaNMDA = deal(rhoCaVDCCslider(1));
%             param.extparam = deal(extslider(1));
             param.Mg_conc = 1;
             param.eta = 0.28;
             param.gamma = 0.016;
            param.NnoiseMit = 0.01*noiseslider(1);
            param.NoiseGradist = 0.01*noiseslider(1);
            param.NoiseGraProx = 0.01*noiseslider(1);
            param.Inoise = 0.1*noiseslider(1);
            
            param.PCinputON = false;
            param.PCtype = pctypeslider(1);
            param.PCparam1 = pcparam1slider(1);
            param.PCparam2 = pcparam2slider(1);

            [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = ...
                        IandVLFP_GCE(Mitral,GraProximal,GraDistal,param);
%%
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

                    scrsz = get(0,'ScreenSize');
                    figH=figure('visible','off'); %not showing all the individual power spectrums
                    set(figH,'position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
                    hold on
                %     plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color',ngradist*cmat(vv,:)); % color
                %     plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color',[ngradist*cmat(vv,:),0.2],'linewidth',2)
                %     plot(f,2*abs(mitFFT1I(1:NFFT/2+1)),'color',ngradist*cmat(vv,:),'linestyle','--')
                %     plot(f,2*abs(mitFFT3I(1:NFFT/2+1)),'color',ngradist*cmat(vv,:),'linestyle',':','linewidth',2);
                    % plot(f,2*abs(mitFFT1V(1:NFFT/2+1)),f,2*abs(mitFFT3V(1:NFFT/2+1)));
                    plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color','k'); % bw
                %   plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color',[0.75,0.75,0.75],'linewidth',2);
                   % plot(f,2*abs(mitFFT1I(1:NFFT/2+1)),'color','k','linestyle','-.');
                   % plot(f,2*abs(mitFFT3I(1:NFFT/2+1)),'color','k','linestyle',':','linewidth',2);

                    %plot(f,2*abs(mitFFTEXTRA(1:NFFT/2+1)),'color','k','linestyle','--'); %[[Sam]]
                    plot(f,2*abs(gradFFTGI(1:NFFT/2+1)),'color','k','linestyle',':'); %[[Sam]]
                    plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color','k','linestyle','-');
                    plot(f,1e-4*ones(1,length(f)),'k--')

                    hold off
                    xlim(faxis(1,:));ylim([0 2e-3])
                    set(gca,'fontsize',fs)
                    %legend('ILFP - Global','VLFP - Global','ILFP - highest 15 MCs','ILFP - lowest 15 MCs')
                    legend('ILFP','ELFP','VLFP') %[Sam]
                    legend boxoff
                    xlabel('Frequency (Hz) ')
                    %title(['PCparam1 ',num2str(param.PCparam1),' Extparam',num2str(param.extparam)])
                    
                    title([slider.name{1},'=',num2str(slider.value{1}(i)),' ',slider.name{2},'=',num2str(slider.value{2}(j))])
                    
                    T.(fieldnames{i}).(slider.name{2})(j) = slider.value{2}(j);
                    [maxpower,ind] = max(2*abs(gradFFTGI(1:NFFT/2+1))); 
                    T.(fieldnames{i}).f(j) = f(ind);
                    T.(fieldnames{i}).power(j) = maxpower;

                    %%
%                     masterdata(i).Mitral = Mitral;
%                     masterdata(i).GraDistal = GraDistal;
%                     masterdata(i).GraProximal = GraProximal;
%                     masterdata(i).InputCurrent = InputCurrent;
%                     masterdata(i).MitLFPs  = MitLFPs;
%                     masterdata(i).GraDistLFPs = GraDistLFPs;
      %  end
    end
end
        %%
%        figure
%        pspectrum(detrend(MitLFPs.GradistMitGlobal(trim:end-100)),10000,'spectrogram','FrequencyLimits',[0 100])
        %%
%         L2 = 1000;
%         NFFT2 = 2^nextpow2(L2);
%         f = sampf/2*linspace(0,1,NFFT2 / 2+1);
%         partialmitFFTGI = fft(detrend(MitLFPs.GradistMitGlobal(2000:3000),'constant'),NFFT2)/L2;
%         figure
%         plot(f,2*abs(partialmitFFTGI(1:NFFT2/2+1)))
%         xlim(faxis(1,:))

%% plot f vs slider and power vs slider
figure('Units','normalized','pos',[0.1,0.1,0.4,0.9])
subplot1 = subplot(2,1,1,'Units','normalized','pos',[0.1,0.52,0.8,0.45]);
subplot2 = subplot(2,1,2,'Units','normalized','pos',[0.1,0.03,0.8,0.45]);
for i=1:length(slider.value{1})
    subT = T.(fieldnames{i});
    axes(subplot1);hold on;plot(subT.(slider.name{2}),subT.f);figname = [slider.name{2},' vs f'];hold off;
    axes(subplot2);hold on;plot(subT.(slider.name{2}),subT.power);figname = [slider.name{2},' vs power'];hold off;
end
axes(subplot1);legend(fieldnames);title([slider.name{2},' vs f']);
axes(subplot2);legend(fieldnames);title([slider.name{2},' vs power']);
savefig([slider.name{1},'_vs_',slider.name{2},'_vs_f_and_power'])
savefile = [slider.name{1},'_vs_',slider.name{2},'_vs.mat'];
save(savefile,'T')
%% turning on datacursor for all figs
figHandles = findobj('type','figure');
for i=1:length(figHandles)
    datacursormode(figHandles(i))
end
