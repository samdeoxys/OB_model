%%%%%%%% Script for parameter sweep. 
%%%%%%%% ----Sam Zheng----

%%%%%%%% To do a single-parameter sweep, create/uncomment the vector of
%%%%%%%% possible values for the parameter in the "Sliders" section, assign
%%%%%%%% the vector to slider.value and its name to slider.name
%%%%%%%% at the end of the "Sliders" section, find/create
%%%%%%%% the relevant parameter in the "Update Parameters" section, change
%%%%%%%% the index on the right hand side from 1 to i, and the current i (if the index for another slider in this secion is i) to
%%%%%%%% 1. This will return the powerspectrum for each value, and 

%%%%%%%% To do a multiple-parameter sweep, create a matrix by stacking multiple sliders with equal
%%%%%%%% length as rows, the rest is similar to single-parameter sweep,
%%%%%%%% except there will be more than one indices to change from 1 to i.
%%%%%%%% This will return the power spectrums but not the max power and
%%%%%%%% frequency plots.

%%%%%%%% A struct array, masterdata, will be returned containing all
%%%%%%%% updated structs for each parameter value (/value pair)

%%%%%%%% If analyze == true, it will return addition figures for analysis.
%%%%%%%% If p_and_f== true, a plot of maximal power vs value and a plot of peak frequency vs value will be returned (only for single-parameter sweep).

%% Global Initialization
clear Mitral GraProximal GraDistal param MitLFPs GraDistLFPs
input_file = 'OB_params_GCE.txt';
[dt,tsim,ntp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);
analyze = true;
p_and_f = true;

%% Sliders
% ampaslider = [0.024];%linspace(0.01,0.17,17); %0.024
% nmdaslider = [0.04]; % 0.02
% nmdatauslider = [4]; % 2 has wide pattern if wnmda=0.04
% nmdatau2slider =75; %75
% ampatauslider = 2;%*linspace(0.1,1.4,20);%2
% ampatau2slider = 3;%*linspace(0.7,3,20);  %3
% %tau2slider = 2*tauslider;
 gabaslider = [0.008,0.0075,0.007]; % 0.008
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
gmodslider = [0.22,0.24,0.26];
%wmod_vrestslider = [Vrestslider;wmodslider];
%gmod_gabaslider = [gmodslider;gabaslider];

ratio = 1; %0.3


% slider that needs to be changed whenever looping on a different variable
slider.value = gmodslider; %Vrestslider; %
slider.name = 'gmod';%[num2str(ratio),' of vrest'];
%% Starting the Loop
masterdata = struct;
T = table();
for i=1:length(slider.value)
            [Mitral GraProximal GraDistal param] = OB_network_GCE(input_file);
            %% Update parameters
            % for changing Vrest for only a subset of GCs
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
               [Mitral.wGABAGR] = deal(gabaslider(1));
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
            %% Run and plot
            [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = ...
                        IandVLFP_GCE(Mitral,GraProximal,GraDistal,param);
            if sum(InputCurrent.CCaBase > GraDistal(1).CCaTh)>0 % warning when the model breaks down
                "CCaBase bigger than CCaTh, something is wrong"
            end
            [figH,mitFFTGV,NFFT,f] = plot_power(MitLFPs,GraDistLFPs,sampf,timevec);
             if size(slider.value,1)==1 % if not it means multiple sweep; only generate T, the table for parameter vs peak power and frequency, when single sweep
                  fig_title = [slider.name,'=',num2str(slider.value(i))]; 
                  T.(slider.name)(i) = slider.value(i);
                  [maxpower,ind] = max(2*abs(mitFFTGV(1:NFFT/2+1))); 
                  T.f(i) = f(ind);
                  T.power(i) = maxpower;
                  title(fig_title)

             else
                 fig_title = [slider.name,'=',mat2str(slider.value(:,i)')];
                 title(fig_title)
             end
            %% storing sweep results in masterdata
            masterdata(i).Mitral = Mitral;
            masterdata(i).GraDistal = GraDistal;
            masterdata(i).GraProximal = GraProximal;
            masterdata(i).InputCurrent = InputCurrent;
            masterdata(i).MitLFPs  = MitLFPs;
            masterdata(i).GraDistLFPs = GraDistLFPs;
            %% analysis
            if analyze == true
                taudist = 100;
                %analyze_firing(Mitral,fig_title); not in use
                analyze_synchrony(Mitral,GraDistal,taudist,fig_title);
                plot_currents_raster(InputCurrent,Mitral,MitLFPs,param.nGradist,timevec,ntp,dt);
            end

end    

%% plot f vs slider and power vs slider
if p_and_f == true && size(slider.value,1)==1
  figure;plot(T.(slider.name),T.f);figname = [slider.name,' vs f.jpg'];title(figname);%saveas(gcf,figname);
  figure;plot(T.(slider.name),T.power);figname = [slider.name,' vs power.jpg'];title(figname);%saveas(gcf,figname);
end
%% turning on datacursor for all figs
figHandles = findobj('type','figure');
for i=1:length(figHandles)
    datacursormode(figHandles(i))
end
allfigname = [slider.name,'_sweep'];
%savefig(figHandles,allfigname)
