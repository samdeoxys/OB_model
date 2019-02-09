clear Mitral GraProximal GraDistal param MitLFPs GraDistLFPs
input_file = 'OB_params_GCE_Fig4.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);
%[Mitral GraProximal GraDistal param] = OB_network_GCE(input_file);


%%
% ampaslider = [0.024];%linspace(0.01,0.17,17); %0.024
% nmdaslider = [0.04]; % 0.02
% nmdatauslider = [4]; % 2 has wide pattern if wnmda=0.04
% nmdatau2slider =75; %75
% ampatauslider = 2;%*linspace(0.1,1.4,20);%2
% ampatau2slider = 3;%*linspace(0.7,3,20);  %3
% %tau2slider = 2*tauslider;
% gabaslider = 0.008; % 0.008
 mittauslider = 7;%linspace(1,10,3);%5; %linspace(1,10,20);
 gcdtauslider = 5;%linspace(1,10,5);
% gcptauslider = 5;%linspace(1,10,20);
% ccathslider = 1.5; %1.5
 extslider = [1.25];%[1.25];
% vdccslider = 200;%linspace(0,400,20); %200
% tauvdccslider = 18;%18*linspace(0.1,0.5,10);%18*linspace(0.25,2,20);%[18];
% collateralslider = [0.2];
% noiseslider = [0];
 taucaslider = [5];%5;
% rhoCaNMDAslider = 50;%50 * linspace(1/2,4,20);
% rhoCaVDCCslider = 100;%50 * linspace(1/2,4,20); % 100

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
etaslider = 0.28;% linspace(0.05,0.25,5);
Vgcdslider = -0.074;
Vgcpslider = -0.074;
Vrestslider = -0.074;%[-0.06,-0.055,-0.05,-0.04,-0.03];%[-0.074 -0.054];
wmodslider = 0;%[0.55,0.58,0.6];%[0.62 0];
%wmod_vrestslider = [Vrestslider;wmodslider];
% slider that needs to be changed whenever looping on a different variable
gmodslider = [0.05,0.1,0.15,0.2];

ratio = 1; %0.3

slider.value = gmodslider; %Vrestslider; %
slider.name = 'gmod';%[num2str(ratio),' of vrest'];
%%
masterdata = struct;
T = table();
for i=1:length(slider.value)
            [Mitral GraProximal GraDistal param] = OB_network_GCE(input_file);
            rng(10)
            ind = find(rand(param.nGradist,1)<=ratio);
            rng('default')
           
%              [GraProximal.tau] =deal(gcptauslider(1));
             [Mitral.tau] = deal(mittauslider(1));
             [GraDistal.tau] = deal(gcdtauslider(1));
%               [GraDistal.wNMDAMI] = deal(nmdaslider(1));
% %             % [GraDistal.wGABAGRSP] = deal(0);
%               [GraDistal.wAMPAMI] = deal(ampaslider(1));
%               [GraProximal.wAMPAMI] = deal(ampaslider(1));
%               [GraProximal.wNMDAMI] = deal(nmdaslider(1));
%              [GraDistal.wVDCCMI] = deal(vdccslider(1));
%              [GraDistal.tauVDCC] = deal(tauvdccslider(1));
%               [Mitral.wGABAGR] = deal(gabaslider(1));
               [GraDistal(ind).Vrest] = deal(Vrestslider(1));
              [GraProximal(ind).Vrest] = deal(Vrestslider(1));
% %              %[GraProximal.FThresh] = deal(-0.05);
%               [GraDistal.tauAMPA1] = deal(ampatauslider(1));
%               [GraProximal.tauPROX1] = deal(1/2*ampatauslider(1));
%               [GraDistal.tauAMPA2] = deal(ampatau2slider(1));
%               [GraProximal.tauPROX2] = deal(2/3*ampatau2slider(1));
%               [GraDistal.tauNMDA1] = deal(nmdatauslider(1));
%               [GraProximal.tauNMDA1] = deal(nmdatauslider(1));
%               [GraDistal.tauNMDA2] = deal(nmdatau2slider(1));
%               [GraProximal.tauNMDA2] = deal(nmdatau2slider(1));
%               param.CChanceMitGraProx = deal(collateralslider(1));
%              [GraDistal.CCaTh] = deal(ccathslider(1)); %1.5
%              [GraDistal.CCaR] = deal(0);
              [GraDistal.tauCa] = deal(taucaslider(1));
% %             
%              param.rhoCaVDCC = deal(rhoCaNMDAslider(1));
%              param.rhoCaNMDA = deal(rhoCaVDCCslider(1));
             param.extparam = deal(extslider(1));
%              param.Mg_conc = 1;
              param.eta = etaslider(1);
%              param.gamma = 0.016;
%             param.NnoiseMit = 0.01*noiseslider(1);
%             param.NoiseGradist = 0.01*noiseslider(1);
%             param.NoiseGraProx = 0.01*noiseslider(1);
%             param.Inoise = 0.1*noiseslider(1);
            
            param.PCinputON = false;
            param.PCtype = pctypeslider(1);
            param.PCparam1 = pcparam1slider(1);
            param.PCparam2 = pcparam2slider(1);
            param.PCparam3 = pcparam3slider(1);
            param.PCparam4 = pcparam4slider(1);
            
            param.wMod = wmodslider(1);
            param.gMod = gmodslider(i);
            
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
                    figH=figure;
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
                     if size(slider.value,1)==1 % if so it means multiple sweep
                           fig_title = [slider.name,'=',num2str(slider.value(i))];
%                          title(num2str(slider.value(i)))
% 
%                          T.(slider.name)(i) = slider.value(i);
%                          [maxpower,ind] = max(2*abs(gradFFTGI(1:NFFT/2+1))); 
%                          T.f(i) = f(ind);
%                          T.power(i) = maxpower;
                           title(fig_title)

                     else
                         fig_title = [slider.name,'=',mat2str(slider.value(:,i)')];
                         title(fig_title)
                     end
                     if sum(InputCurrent.CCaBase > GraDistal(1).CCaTh)>0
                         "CCaBase bigger than CCaTh, something is wrong"
                     end
                    
                    %%
                    masterdata(i).Mitral = Mitral;
                    masterdata(i).GraDistal = GraDistal;
                    masterdata(i).GraProximal = GraProximal;
                    masterdata(i).InputCurrent = InputCurrent;
                    masterdata(i).MitLFPs  = MitLFPs;
                    masterdata(i).GraDistLFPs = GraDistLFPs;
      %  end
  %  end
  
  %% analysis
                    taudist = 100;
                    analyze_firing(Mitral,fig_title);
                    analyze_synchrony(Mitral,GraDistal,taudist,fig_title);
                    plot_currents_raster(InputCurrent,Mitral,MitLFPs,param.nGradist,timevec,ntp,dt);
end
    

%% plot f vs slider and power vs slider
 % figure;plot(T.(slider.name),T.f);figname = [slider.name,' vs f.jpg'];title(figname);saveas(gcf,figname);
 % figure;plot(T.(slider.name),T.power);figname = [slider.name,' vs power.jpg'];title(figname);saveas(gcf,figname);
%% turning on datacursor for all figs
figHandles = findobj('type','figure');
for i=1:length(figHandles)
    datacursormode(figHandles(i))
end
allfigname = [slider.name,'_sweep'];
%savefig(figHandles,allfigname)
