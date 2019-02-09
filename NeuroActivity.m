%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mitral GraProximal GraDistal param InputCurrent] = NeuroActivity(Mitral,GraProximal,GraDistal,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates the neuronal activity
%
% Licurgo de Almeida
% 11/03/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create initial setup -- OB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

if param.flagRespiration == true
    Respiration = CreateRespFreq(param);
else
    Respiration = ones(1,round(param.tsim / param.dt));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%   OB   %%%%%

nds = param.nGradist/param.nGraprox; % number of distal synapses per granule cell
%????


    % Each MC connected to random subset CChanceGraMit*nGradist dGCs
    % In this scheme each MC is connected to exactly CChanceGraMit*nGradist
    % dGCs, but each dGC can be connected to a range of MCs
        rng(2)
        MatGradistMit = SetConnections(param.nMitral,param.nGradist,param.CChanceGraMit);
        
        MatMitGradist = MatGradistMit'; % reciprocal synapses
        
        rng(3)
        MatMitGraprox = SetConnections(param.nGraprox,param.nMitral,param.CChanceMitGraProx); %[[Sam]] connection mat for MC to graprox, for axon collateral
        MatGraproxMit = MatMitGraprox';
        rng(4) %%[[Sam]] to be removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set weights matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%   OB   %%%%%

    wGABAMit = zeros(param.nMitral,1);
    wGABASPMit = zeros(param.nMitral,1);
    for ii = 1:param.nMitral
        wGABAMit(ii) = Mitral(ii).wGABAGR;
        wGABASPMit(ii) = Mitral(ii).wGABAGRSP; %%[[Sam]] for gc spike dependent gaba
    end
    wGradistMit = SetWeights(MatGradistMit,wGABAMit);
    
    wGraproxMit = SetWeights(MatGraproxMit,wGABASPMit); %%[[Sam]] for gc spike dependent gaba
    
    wAMPAGradist = zeros(param.nGradist,1);
    wNMDAGradist = zeros(param.nGradist,1);
    wVDCCGradist = zeros(param.nGradist,1);
    wGABAGraProx = zeros(param.nGraprox,1);
    wIntraGradist = zeros(param.nGradist,1); %[Sam]
    wIntraGraProx = zeros(param.nGraprox,1); %
    for ii = 1:param.nGradist
        wAMPAGradist(ii) = GraDistal(ii).wAMPAMI;
        wNMDAGradist(ii) = GraDistal(ii).wNMDAMI;
        wVDCCGradist(ii) = GraDistal(ii).wVDCCMI;
        wIntraGradist(ii) = GraDistal(ii).wIntraDist; %[Sam]
        
    end
    for ii = 1:param.nGraprox
        wGABAGraProx(ii) = GraProximal(ii).wGABAGR;
        wAMPAGraProx(ii) = GraProximal(ii).wAMPAMI; % [[Sam]]
        wIntraGraProx(ii) = GraProximal(ii).wIntraProx; % [[Sam]]
        wNMDAGraProx(ii) = GraProximal(ii).wNMDAMI; %[[Sam]] added nmda
    end
    wMitGradistAMPA = SetWeights(MatMitGradist,wAMPAGradist);
    wMitGradistNMDA = SetWeights(MatMitGradist,wNMDAGradist);
    wMitGradistVDCC = SetWeights(MatMitGradist,wVDCCGradist);
    %MatMitGraProx = MatMitGradist;
    wMitGraProxAMPA = SetWeights(MatMitGraprox,wAMPAGraProx); % [[Sam]]
    wMitGraProxNMDA = SetWeights(MatMitGraprox,wNMDAGraProx); % [[Sam]]
%     wProxProx = SetWeights(MatProxProx,wGABAGraProx);
    
    if strcmp(param.PCtype,'poisson')
        wpcinput = param.PCparam2*wAMPAGraProx .* eye(param.nGraprox); %[7.0] weight of poisson pcinput; important: .*!; diagnal matrix so that in SetI dimensions are correct
    else
        %wpcinput = param.PCparam2 * eye(param.nGraprox).*sort(rand(param.nGraprox,1),'descend'); %[7.0] the problem with weight coded with uniform distribution is the largest and smallest really different; consider gaussian with different mean later!
        wpcinput = param.PCparam2 * eye(param.nGraprox); % [Sam]without randomness
    end
    
    %[7.0] PCparam2 in the poisson case is for scaling
    %param.wpcinput = wpcinput; % for external access and debugging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Mitral cell at a given time
Vmit = zeros(param.nMitral,round(param.tsim / param.dt));
Smit = Vmit; % Binary matrix recording only the spikes
Iext_matrix = Vmit; % Initialize external inputs to MCs, (all 0 for now)

refracmit = 3; % Refractory period after firing (ms)

Restmit = zeros(param.nMitral,1); % resting potential   ???why 0 and why is it an array?
Threshmit = Restmit; % current firing threshold
Hypermit = Restmit; % hyperpolarization potential  ??? why all Restmit??? to initialize?? not necessary?
taumit = Restmit; % tau neuron
countrefracmit = Restmit; % Refractory period counter. This is assumed to
% be the same for all Mitral cells so we don't need to add it to the for
% below.
gmaxAMPAmit = Restmit; % Max AMPA conductance
tauAMPA1mit = Restmit; % AMPA's rising tau.
tauAMPA2mit = Restmit; % AMPA's falling tau.
EAMPAmit = Restmit; % AMPA's Nernst potential.
gmaxGABAmit = Restmit; % Max GABA conductance
tauGABA1mit = Restmit; % GABA's rising tau.
tauGABA2mit = Restmit; % GABA's falling tau.
EGABAmit = Restmit; % GABA's Nernst potential.


for ii = 1:param.nMitral
    Restmit(ii) = Mitral(ii).Vrest;
    Hypermit(ii) = Mitral(ii).Vhyper;
    taumit(ii) = Mitral(ii).tau;
    gmaxAMPAmit(ii) = Mitral(ii).gmaxAMPA;
    tauAMPA1mit(ii) = Mitral(ii).tauAMPA1;
    tauAMPA2mit(ii) = Mitral(ii).tauAMPA2;
    EAMPAmit(ii) = Mitral(ii).EAMPA;
    gmaxGABAmit(ii) = Mitral(ii).gmaxGABA;
    tauGABA1mit(ii) = Mitral(ii).tauGABA1;
    tauGABA2mit(ii) = Mitral(ii).tauGABA2;
    EGABAmit(ii) = Mitral(ii).EGABA;

    Threshmit(ii) = Mitral(ii).FThresh; 
    
    Mitral(ii).ConnectionsGrad = MatGradistMit(ii,:);
    Mitral(ii).ConWeightsGrad = wGradistMit(ii,:);
    
    Mitral(ii).ConnectionsGrap =MatGraproxMit(ii,:);
    Mitral(ii).ConWeightsGrad  = wGraproxMit(ii,:);
end

% Initialize Mitral cells potentials
Vmit(:,1) = Restmit;
Vmit_nospike = Vmit;


% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Only used if paramal.ProxON is
% true and param.DistalON is false (i.e. spiking inhibition)
tGABA0mit = zeros(param.nGradist,round(param.tsim / param.dt)) - 10000000; %?????

Igradistmit = zeros(param.nMitral,1); % Input coming from Granule cells
Igradistmit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));

Vnoisemit = zeros(param.nMitral,1); %[[Sam]]test when no noise
Vnoisemit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));

% only used when distal gra dendrites are removed
Igraproxmit = zeros(param.nMitral,1); % Input coming from Granule cells
Igraproxmit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));

% [[Sam]] for spike dependent gaba
Igraspmit = zeros(param.nMitral,1);
Igraspmit_matrix = zeros(param.nMitral,round(param.tsim/param.dt));

maxgGABAmit = getmaxg(param.dt,tauGABA1mit,tauGABA2mit); % Get max conductance
% amplitude


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set Proximal Granule cell parameters and variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.ProximalON == true
  Vgraprox = zeros(param.nGraprox,round(param.tsim / param.dt));
  Sgraprox = Vgraprox;
% % Probability of each cell firing
  Pfiregraprox = Vgraprox;
% 
  refracgraprox = 5; % Refractory period after firing (ms)
% 
  Restgraprox = zeros(param.nGraprox,1); % resting potential
  Threshgraprox = Restgraprox; % current firing threshold
  Hypergraprox = Restgraprox; % hyperpolarization potential
  taugraprox = Restgraprox; % tau neuron
  countrefracgraprox = Restgraprox; % Refractory period counter. This is assumed to
% % be the same for all Granule cells so we don't need to add it to the for below.
% 
  spikelag = Restgraprox; %[[Sam]]
  countlag = Restgraprox; %[[Sam]]
  gmaxAMPAgraprox = Restgraprox; % Max AMPA conductance
% 
  tauPROX1 = Restgraprox; % AMPA's rising tau.
  tauPROX2 = Restgraprox; % AMPA's falling tau.
% 
  EAMPAgraprox = Restgraprox; % AMPA's Nernst potential.
  gmaxGABAgraprox = Restgraprox; % Max GABA conductance
  tauGABA1graprox = Restgraprox; % GABA's rising tau.
  tauGABA2graprox = Restgraprox; % GABA's falling tau.
  EGABAgraprox = Restgraprox; % GABA's Nernst potential.
  ENMDAgraprox = Restgraprox; % [[Sam]]
% 
% 
  for ii = 1:param.nGraprox
      Restgraprox(ii) = GraProximal(ii).Vrest;
      Hypergraprox(ii) = GraProximal(ii).Vhyper;
      taugraprox(ii) = GraProximal(ii).tau;
      gmaxAMPAgraprox(ii) = GraProximal(ii).gmaxAMPA;
      tauPROX1(ii) = GraProximal(ii).tauPROX1;
      tauPROX2(ii) = GraProximal(ii).tauPROX2;
      tauNMDA1(ii) = GraProximal(ii).tauNMDA1; %[[Sam]] added nmda
      tauNMDA2(ii) = GraProximal(ii).tauNMDA2;
      EAMPAgraprox(ii) = GraProximal(ii).EAMPA;
      gmaxGABAgraprox(ii) = GraProximal(ii).gmaxGABA;
      tauGABA1graprox(ii) = GraProximal(ii).tauGABA1;
      tauGABA2graprox(ii) = GraProximal(ii).tauGABA2;
      EGABAgraprox(ii) = GraProximal(ii).EGABA;
      ENMDAgraprox(ii) = GraProximal(ii).EAMPA;
      
      Threshgraprox(ii) = GraProximal(ii).FThresh; %[[Sam]] proximal firing threshold, for spike dependent gaba
    %  Threshgraprox(ii) = GraProximal(ii).CCaTh; %[Sam]does proximal has any
    %  use for CCaTh?
      GraProximal(ii).Connectionsmit =MatGraproxMit(:,ii)';
      spikelag(ii) = GraProximal(ii).spikelag; %[Sam] the time it has to be above threshold to fire;
  end
% 
% % Initialize Granule cells potentials
  Vgraprox(:,1) = Restgraprox;
  Vgraprox_nospike = Vgraprox;
% 
% % bulbar input current to graprox
  Idistprox = zeros(param.nGraprox,1); % Input coming from dist to prox [add semicolon]
  Idistprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
  ImitgraproxAMPA = zeros(param.nGraprox,1); %[[Sam]]
  ImitgraproxAMPA_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
  
  ImitgraproxNMDA = zeros(param.nGraprox,1); %[Sam]added nmda
  ImitgraproxNMDA_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
  
  Mg_block_graprox = zeros(param.nGraprox,1); %[Sam]
  Mg_blcok_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
  
  Vnoisegraprox = zeros(param.nGraprox,1); %[[Sam]]
  Vnoisegraprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
% 
% % input voltage (for direct distprox summing)
  Vdistprox = zeros(param.nGraprox,1); % V input coming from dist to prox
  Vdistprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
% 
% 
% % AMPA time counter. This variable starts with a very negative value just
% % to make sure that the currents will be = 0
% 
% % tAMPA0mit_gra is a matrix that stores tau values delayed by mit-gra delay time
  tAMPA0mit_gra = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;
% 
% 
  maxgAMPAgraprox = getmaxg(param.dt,tauPROX1,tauPROX2); % Get max conductance amplitude
  maxgNMDAgraprox = getmaxg(param.dt,tauNMDA1,tauNMDA2); % Get max conductance amplitude [[Sam]]
% 
% % maxgAMPAgraprox = 1; % Set gmax = 1 because it will be normalized anyways
% 
% % GABA time counter. This variable starts with a very negative value just
% % to make sure that the currents will be = 0. Each Granule cell can be
% % connected to a varied number of other granule cells
 tGABA0graprox = zeros(param.nGraprox,round(param.tsim / param.dt)) - 10000000;
 Igragra = zeros(param.nGraprox,1); % Input coming from other granule cells
 Igragra_matrix = zeros(param.nGraprox,round(param.tsim / param.dt)); %uncommented so that it is defined in calculating graprox membrane potential
% 
% maxgGABAgraprox = getmaxg(param.dt,tauGABA1graprox,tauGABA2graprox); % Get max conductance amplitude
% 
% % Sam: PC input to graprox
  Ipcgraprox = zeros(param.nGraprox,1);
  Ipcgraprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
  if param.PCinputON == true && strcmp(param.PCtype,'poisson') == true %[7.0] adding PCspike for poisson input
      pcspikes = zeros(param.nGraprox);
      pcspikes_matrix = zeros(param.nGraprox,round(param.tsim/param.dt));

  end
t0pc_gra = zeros(param.nGraprox,round(param.tsim/param.dt)) - 10000000; % storing time of last firing; had to initialize even if not poisson, so that it can be passed to the function SetI_PC
      
  
  
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % End of Proximal Granule cells parameters and variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Distal Granule cell parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.DistalON == true
% Stores the voltage of each Granule cell at a given time
Vgradist = zeros(param.nGradist,round(param.tsim / param.dt));
CCa = zeros(param.nGradist,round(param.tsim / param.dt));% Ca concentration
taum_matrix = CCa;
Sgradist= Vgradist;% Binary matrix recording only the spikes


Prelease_matrix = Vgradist; % Prelease_matrix is the probability of GABA release 
                             % when using graded inhibition from distal gra dendrites

Restgradist = zeros(param.nGradist,1); % resting potential
Threshgradist = Restgradist; % current firing threshold
Hypergradist = Restgradist; % hyperpolarization potential
taugradist = Restgradist; % tau neuron
gmaxAMPAgradist = Restgradist; % Max AMPA conductance


tauAMPA1 = Restgradist; % AMPA's rising tau.
tauAMPA2 = Restgradist; % AMPA's falling tau.
tauNMDA1 = Restgradist; % NMDA's rising tau.
tauNMDA2 = Restgradist; % NMDA's falling tau.
tauVDCC = Restgradist; % VDCC activation tau.

tauCagradist = Restgradist; % decay of ca

% L-type Ca activation parameters
mbar = Restgradist;
taum = Restgradist;
m = Vgradist;
h = Vgradist;

EAMPAgradist = Restgradist; % AMPA's Nernst potential.
ENMDAgradist = Restgradist; % NMDA's Nernst potential.
EVDCCgradist = Restgradist; % Ca Nernst potential.
CCaR = Restgradist; % initialize to make the dimensions right!!

for ii = 1:param.nGradist
    Restgradist(ii) = GraDistal(ii).Vrest;
    Hypergradist(ii) = GraDistal(ii).Vhyper;
    taugradist(ii) = GraDistal(ii).tau;
    gmaxAMPAgradist(ii) = GraDistal(ii).gmaxAMPA;
    tauAMPA1(ii) = GraDistal(ii).tauAMPA1;
    tauAMPA2(ii) = GraDistal(ii).tauAMPA2;
    tauNMDA1(ii) = GraDistal(ii).tauNMDA1;
    tauNMDA2(ii) = GraDistal(ii).tauNMDA2;
    tauVDCC(ii) = GraDistal(ii).tauVDCC;
    EAMPAgradist(ii) = GraDistal(ii).EAMPA;
    ENMDAgradist(ii) = GraDistal(ii).EAMPA; % nmda reversal potential is the same as AMPA
    EVDCCgradist(ii) = GraDistal(ii).ECa; % Ca reversal potential is ~100mV higher than AMPA
    % EVDCCgradist(ii) = GraDistal(ii).EAMPA;
    
    Threshgradist(ii) = GraDistal(ii).CCaTh;
    tauCagradist(ii) = GraDistal(ii).tauCa; %[Sam] tau of Ca decay
    CCaR(ii) = GraDistal(ii).CCaR; %[Sam] 
    
    GraDistal(ii).Connections = MatMitGradist(ii,:);
    
end
% 
% 1-ExFrac subpopulation of GCs in unexcited state
if param.ExFrac < 1
    Restgradist((round(param.ExFrac*param.nGradist)+1):end) = -75e-3;
end %???

% Initialize Granule cells potentials
Vgradist(:,1) = Restgradist;
%CCa(:,1) = 0.1; % [uM] ????

% proximal granule input currents to gradist
Iproxdist = zeros(param.nGradist,1); % Input coming from prox to dist
Iproxdist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% mitral input to gradist AMPA and NMDA receptors
ImitgradistAMPA = zeros(param.nGradist,1);
ImitgradistAMPA_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

ImitgradistNMDA = zeros(param.nGradist,1);
ImitgradistNMDA_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

Mg_block_gradist = zeros(param.nGradist,1); %[Sam]
Mg_block_gradist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));
% L-type Ca current
ImitgradistVDCC = zeros(param.nGradist,1);
ImitgradistVDCC_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% Vgradist noise
Vnoisegradist = zeros(param.nGradist,1); %[[Sam]]no noise condition
Vnoisegradist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% Imodgradist
Imodgradist = zeros(param.nGradist,1);
Imodgradist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));


% tAMPA0mit_gra is a matrix that stores tau values delayed by mit-gra delay time
tAMPA0mit_gra = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000; 
% Wmitgradist = ones(param.nGranule,param.nMitral); % Synaptic weights

maxgAMPAgradist = getmaxg(param.dt,tauAMPA1,tauAMPA2); % Get max conductance amplitude for AMPA
maxgNMDAgradist = max(exp(-([0:.01:200])./tauNMDA2(1)) - exp(-([0:.01:200])./tauNMDA1(1))); % Get max conductance amplitude for NMDA


% Calculate baseline [Ca] from steady state solution of [Ca]' = 0
mbarrest = 1./(1 + exp(-(1e3*Restgradist+45)./7)); % Restgradist must be in units of mV
RT = 300*8.31; % J/mole
z= 2 ; % Ca ion valence
F = 96485; % Faraday constant Coul/mole
Cout = 1500;
Cin = 0.002:0.002:2.5; %???
EVDCCgradist = (RT/(z*F))*log(Cout./Cin); % RHS of equation
LHS = zeros(param.nGradist,length(Cin));
for ii = 1:length(Cin)
    for nn = 1:param.nGradist
        LHS(nn,ii) = Restgradist(nn) + (Cin(ii)^2)./(1e-4 .* param.rhoCaVDCC .* sum(wMitGradistVDCC(nn,:)) .* mbarrest(nn)); %??? %denominator could be zero in certain circumstances, careful; 
        
    end
end

% plot LHS and RHS
% subplot(2,1,1)
% plot(Cin,LHS(1,:),'k.',Cin,EVDCCgradist,'k--')
% set(gca,'fontsize',16)
% legend('LHS','RHS','location','best')
% % plot |LHS - RHS|
% subplot(2,1,2)
% plot(Cin,abs(LHS(1,:) - EVDCCgradist),'k')
% set(gca,'fontsize',16)
% title('abs(LHS - RHS)')
% xlabel('[Ca] (\muM)')

% old solution (without explicit [Ca] dpendence in Eca)
% CCaBase = sqrt(1e-4 * nmit_per_gra .* param.rhoCa .* wVDCCGradist .* mbarrest .* (EVDCCgradist - Restgradist));
% Note: the 1e-4 comes from the definition of h for N-type current:
%    h = 1e-4/(1e-4 + [Ca])

% new solution (with explicit [Ca] dpendence in Eca)
CCaBase = zeros(param.nGradist,1);
if wVDCCGradist > 0
    for ii = 1:param.nGradist
    % mind = find(abs(LHS(ii,:) - EVDCCgradist) == min(abs(LHS(ii,:) - EVDCCgradist)));

    
    index = find(abs(LHS(ii,:) - EVDCCgradist) == min(abs(LHS(ii,:) - EVDCCgradist))); %[[Sam]] took what's previously in Cin(), and assigned it to index; doesnt really matter
    if sum(wMitGradistVDCC(ii,:))~=0 %prevent bug caused by gc with no connected mc
        try
        CCaBase(ii) = Cin(index);
        catch
            ii
            abs(LHS(ii,:) - EVDCCgradist)
            min(abs(LHS(ii,:) - EVDCCgradist))
            find(abs(LHS(ii,:) - EVDCCgradist) == min(abs(LHS(ii,:) - EVDCCgradist)))
        end
    else
        CCaBase(ii) = 0;
    end
end
end


CCa(:,1) = 0.1; %[[Sam]]initializing CCaBase

%%%%%%%%????%%%%%%%%%
if sum(wVDCCGradist) > 0
    EVDCCgradistBase = (RT/(z*F))*log(Cout./CCaBase);
    WMGVDCCsum = zeros(param.nGradist,1);
    for ii = 1:param.nGradist
        WMGVDCCsum(ii) = sum(wMitGradistVDCC(ii,:));
    end
    IVDCCBase = 1e-4*WMGVDCCsum.*mbarrest.*(EVDCCgradistBase - Restgradist)./(1e-4 + CCaBase);
else
    IVDCCBase = zeros(param.nGradist,1);
end

end

%[[Sam]]helper variables for numerical integration


tD = 0;
tP = 0;
ID = zeros(param.nGradist,1);
IP = zeros(param.nGraprox,1);
wP = 1;
wD = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Distal Granule cell synapse parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin neuron simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iext_max = SetExtInput(param); % set level of external input
% sc = [rand(1,1000).*tanh(linspace(0,pi,1000)) ones(1,round(param.tsim / param.dt)-1000)];
% scaling function designed to soften initial impact


Iext_matrix = Iext_max; %[[Sam]] no need to assign the matrix at each time step

% Start loop
for tt = 2:round(param.tsim / param.dt)
    
    t = tt * param.dt; % current time
    
    

    % Mitral Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Get Glo input to Mitral cells. This current is applied directly in
    % Vmit, no channel conductance.

    % External inupt to Mit cells is modulated by respiration if Resp is on
%    Iext = Respiration(tt) .* Iext_max(:,tt);
    
    Iext = Iext_matrix(:,tt); %[[Sam]] Included respiration effect in the function for setting external I
%     Iext = sc(tt) * Respiration(tt) .* Iext_max(:,tt);
 %   Iext_matrix(:,tt) = Iext(:);
    
            
    % Get Granule graded inputs to Mitral cells
    Igradistmit(:) = SetInoSpike_GraMit(gmaxGABAmit,Prelease_matrix(:,tt-1),EGABAmit,Vmit_nospike(:,tt - 1),...
        wGradistMit);
    Igradistmit_matrix(:,tt) = Igradistmit(:); 
    Vnoisemit = param.noisemit .* randn(param.nMitral,1);
    Vnoisemit_matrix(:,tt) = Vnoisemit(:);
    
    
    % [[Sam]] for spike dependent gaba to MC
    Igraspmit(:) = SetI(tauGABA1mit(1),tauGABA2mit(1),t,tGABA0mit(:,tt),maxgGABAmit(1),...
       wGradistMit,EGABAmit(1),Vmit_nospike(:,tt - 1)); % the spike dependent gaba to MC still uses wGradistMit instead of wGraproxMit, which is only used for axon collateral
    Igraspmit_matrix(:,tt) = Igraspmit(:);
    % Mitral cell potential
    
    % Forwards Euler
    
    Vmit(:,tt) = Vmit(:,tt - 1) + (param.dt ./ taumit(:)) .* ...
        ((Iext(:) + Igradistmit(:) + Igraspmit(:)) - Vmit(:,tt - 1) + Restmit(:) + Vnoisemit); 
   
    % Backwards Euler
    
%     g_gramit = (((tauGABA1mit(1) * tauGABA2mit(1)) / (tauGABA1mit(1) - tauGABA2mit(1))) *...
%     (exp(-(t - tGABA0mit(:,tt)) / tauGABA1mit(1)) - exp(-(t - tGABA0mit(:,tt)) / tauGABA2mit(1)))) / maxgGABAmit(1);
% 
%     Vmit(:,tt) = (Vmit(:,tt - 1) + (param.dt ./ taumit(1)) .* ((Iext(:) + (wGraMit_scaled * g_gramit) .* EGABAmit(1) + Restmit(:))))...
%         ./ (1 + (param.dt ./ taumit(1)) .* ((wGraMit_scaled * g_gramit) + 1));
   
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vmit(:,tt - 1) == param.SpikeV; %pick out the ones that spike at tt-1; so impossible for a step that's too big? (tt-2 smaller than threshold, tt-1 bigger;)
    Vmit(I,tt) = Hypermit(I); %?? is hyperpolerization potential a fixed value for a particular neuron? so that if it does not spike, it is set to that value
    Vmit_nospike(:,tt) = Vmit(:,tt); %??what does _nospike mean?
    
    % I is a vector of 1s or 0s of length ncells
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Granule variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This variable is used for input onto distal gra dendrites
        Mit_spike_time = t; %%real time, instead of tt; 
        tAMPA0mit_gra(I,ceil(Mit_spike_time/param.dt):end) = Mit_spike_time; %%??? 
        tNMDA0mit_gra = tAMPA0mit_gra;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracmit(I) = refracmit / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracmit > 0; %%so now I is not just those that fire at tt-1 but also anyone still in refractory period
    countrefracmit(I) = countrefracmit(I) - 1; %%update the refractory time left;
    
    I = find(countrefracmit == 0); % if countrefracmit = 0 the neuron can fire
    % again

% spike when V > thresh    
    if ~isempty(I)
    J = find(Vmit(I,tt) >= Threshmit(I));
    if ~isempty(J)
            Vmit(I(J),tt) = param.SpikeV; % Action potential
            Smit(I(J),tt) = 1; % Record spike time
    end
    end
    
    
    % Granule Distal Dendrites
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% [[Sam_assumption]] now assume distal always on, and prox can only
    %%% be on when distal is on; so proximal condition is nested within the
    %%% distal condition
if param.DistalON == true
       
    % NOTE: In the paper the weights were not included in the definition of the 
    % currents. But we include them in the calculation of the currents because
    % multiplication by weight matrix allows convenient summation of all
    % inputs onto cell. The result is the same.
    
    
    %%%% currents and other variables for gradistal %%%%
    wD = GraDistal(1).wIntraDist; %[[Sam]]helper variables for integration
    tD = param.dt/taugradist(1);
    
    % Get Mitral input to Distal Granule AMPA receptors
    ImitgradistAMPA(:) = SetI(tauAMPA1(1),tauAMPA2(1),t,tAMPA0mit_gra(:,tt),maxgAMPAgradist(1),...
        wMitGradistAMPA,EAMPAgradist(1),Vgradist(:,tt - 1));
    ImitgradistAMPA_matrix(:,tt) = ImitgradistAMPA(:);
    
    % Mitral input to Distal Granule NMDA receptors
    [ImitgradistNMDA(:), Mg_block_gradist(:)] = SetI_NMDA(tauNMDA1,tauNMDA2,t,tNMDA0mit_gra(:,tt),maxgNMDAgradist,...
        wMitGradistNMDA,ENMDAgradist(1),Vgradist(:,tt - 1),param);
    ImitgradistNMDA_matrix(:,tt) = ImitgradistNMDA(:);
    Mg_block_gradist_matrix(:,tt) = Mg_block_gradist(:);
    % Distal Granule VDCC Current
    
            % Calcium Concentration variable
    CCa(:,tt) = CCa(:,tt - 1) + (param.dt ./ tauCagradist(:)) .* (param.rhoCaNMDA .*...
          ImitgradistNMDA(:) + param.rhoCaVDCC.*ImitgradistVDCC(:) ...
          - CCa(:,tt - 1) );
      
%     % L-type
%     % % m' = (mbar - m)/taum;
%    mbar(:) = 1./(1 + exp(-(1e3.*Vgradist(:,tt - 1)+50)./3));
%    taum(:) = 18 * exp(-((1e3.*Vgradist(:,tt - 1)+45)./20).^2) + 1.5;
%    h(:,tt) = 0.00045./(0.00045+CCa(:,tt));

    % Ntype
    mbar(:) = 1./(1 + exp(-(1e3.*Vgradist(:,tt - 1)+45)./7)); %% [Sam]changing the half activation value from 45!
    taum(:) = tauVDCC .* exp(-((1e3.*Vgradist(:,tt - 1)+70)./25).^2) + 0.3; %% [Sam] changing denominator term from 25
    if param.hCaflag == true
     h(:,tt) = 0.0001./(0.0001+CCa(:,tt));
    else
        h(:,tt) = 0.0001./(0.0001+CCaBase);
        % h(:,tt) = 0.0001./(0.0001+WMGVDCCsum*2.0714e-04);
        % scaling factor to obtain mean [Ca] for WVDCC = 350, WGM = 0.025
        % is 2.0714e-04
    end
     
    mtemp = m(:,tt - 1) + (param.dt ./ taum(:)).*(mbar(:) - m(:,tt - 1));
    mtemp(mtemp < 0) = 0;
    m(:,tt) = mtemp;
    taum_matrix(:,tt) = taum(:);
    
%     RT = 300*8.31; % J/mole
%     z= 2 ; % Ca ion valence
%     F = 96485; % Faraday constant Coul/mole
%     Cout = 1500;
    Cin = CCa(:,tt);

    EVDCCgradist = (RT/(z*F))*log(Cout./Cin);

    ImitgradistVDCC(:) = SetI_VDCC(wMitGradistVDCC,EVDCCgradist,Vgradist(:,tt - 1),m(:,tt),h(:,tt));
%     ImitgradistVDCC(ImitgradistVDCC<0) = 0;
    ImitgradistVDCC_matrix(:,tt) = ImitgradistVDCC(:);
    
    
    % NOTE!!! All voltages are in V, not mV, so parameters have
    % to be scaled accordingly!!!

    


    Vnoisegradist = param.noisegradist .* randn(param.nGradist,1);
    Vnoisegradist_matrix(:,tt) = Vnoisegradist;
    
    %
    Imodgradist = param.gMod.*(param.EMod - Vgradist(:,tt-1));
    Imodgradist_matrix(:,tt) = Imodgradist;
   
    %fill in helper variables
    ID(:) = ImitgradistAMPA(:) + ImitgradistNMDA(:) + (ImitgradistVDCC(:) - IVDCCBase) + (1-param.wMod)*(Restgradist(:) +Imodgradist - Vgradist(:,tt - 1)) + Vnoisegradist; %synaptic + leaky + noise current to Gradistal; neuromodulation
    %%%% currents and other variables for gradistal end %%%%

     % Get Proximal Granule input to Distal Granule Dendrites [[[====================Sam==================]]]
     if param.ProximalON == true
         %%%% variables and update V %%%%
         wP = GraProximal(1).wIntraProx; %[Sam] helper variables for numerical integration
         tP = param.dt/taugraprox(1);
         
         %%% turn on/off Imitgraprox, axon collateral
         if param.CollateralON == true
          ImitgraproxAMPA(:) = SetI(tauPROX1(1),tauPROX2(1),t,tAMPA0mit_gra(:,tt),maxgAMPAgraprox(1),...
              wMitGraProxAMPA,EAMPAgraprox(1),Vgraprox_nospike(:,tt-1));
          [ImitgraproxNMDA(:),Mg_block_graprox(:)]= SetI_NMDA(tauNMDA1,tauNMDA2,t,tNMDA0mit_gra(:,tt),maxgNMDAgraprox(1),...
        wMitGraProxNMDA,ENMDAgraprox(1),Vgraprox_nospike(:,tt - 1),param); %[[Sam]] added nmda
         end
         ImitgraproxAMPA_matrix(:,tt) = ImitgraproxAMPA(:);
         ImitgraproxNMDA_matrix(:,tt) = ImitgraproxNMDA(:);
         Mg_block_graprox_matrix(:,tt) = Mg_block_graprox(:);
         
         if param.PCinputON == true && param.ProximalON == true % set PC input to graprox 
            if strcmp(param.PCtype,'poisson') == true %[7.0] adding PCspike for poisson input
                pcspikes = rand(param.nGraprox,1) < param.PCparam1 * (param.dt/1000); %[7.0] for poisson, PCparam1 is the firing rate, times param.dt/1000 to get Pfiring within the interval;
                %pcspikes = rand < param.PCparam1 * (param.dt/1000);pcspikes = repmat(pcspikes,param.nGraprox,1); %second way, synchronous poisson arrival
                I = pcspikes == 1;
                t0pc_gra(I,tt:end) = t;
                pcspikes_matrix(:,tt) = pcspikes(:);
            
                
                Ipcgraprox = SetI(param.PCparam3,param.PCparam4,t,t0pc_gra(:,tt),maxgAMPAgraprox(1),... % poisson I similar same as AMPA prox
                wpcinput,EAMPAgraprox(1),Vgraprox_nospike(:,tt-1)); %now PCparam3 and PCparam4 for rise and fall time
                %[Ipcgraprox,~] = SetI_NMDA(tauNMDA1,tauNMDA2,t,t0pc_gra(:,tt),maxgNMDAgraprox(1),... % poisson try using NMDA prox
                %wpcinput,EAMPAgraprox(1),Vgraprox_nospike(:,tt-1),param);
                
                
            elseif strcmp(param.PCtype,'constant') == true
                Ipcgraprox = wpcinput * (param.PCparam1 * ones(param.nGraprox,1)); %[7.0] weight diagonal matrix dot column vector of strength
                
                %Ipcgraprox = SetI_PC(param,tt); %now have arguments param and tt, and only for one time step
            elseif strcmp(param.PCtype,'sine') == true
                Ipcgraprox = wpcinput * ones(param.nGraprox,1) * (1+param.PCparam1 * sin(2*pi*param.PCparam3/(1000/param.dt)*tt)); %[Sam] PCparam2 is weight coded in wpcinput, PCparam1 is amplitude for the fluctuation, PCparam3 is frequnecy;
            
            end
         end
         Ipcgraprox_matrix(:,tt)=Ipcgraprox(:); % reversed the equality[7.0]
         %Vnoisegraprox = param.noisegraprox .* randn(param.nGraprox,1);
         Vnoisegraprox = 0; %[[Sam]] no noise condition
         Vnoisegraprox_matrix(:,tt) = Vnoisegraprox(:);
         IP(:) = ImitgraproxAMPA(:) + Ipcgraprox(:) + (1-param.wMod)*(Restgraprox(:) - Vgraprox(:,tt-1)) + Vnoisegraprox + ImitgraproxNMDA(:);
         
         Vgraprox(:,tt) = tP*wP/(1+tP*wP+tD*wD) .* (Vgradist(:,tt-1) + tD.*ID(:) + ((1+tD*wD)/wP).*IP(:) + ((1+tD*wD)/(tP*wP)*Vgraprox(:,tt-1)));
         %%%% variables and update V end %%%%
         
         
         %%%% post V change effects, spike %%%%
         % If the neuron fired last cycle, neuron potential hyperpotentializes
        
         I = Vgraprox(:,tt-1)==param.SpikeV;
         Vgraprox(I,tt) = Hypergraprox(I);
         Vgraprox_nospike(:,tt) = Vgraprox(:,tt);
         %%%%%%%%%%%%%%%%[[Sam]] for GC spike dependent gaba onto MC
            GraProx_spike_time = t;
            tGABA0mit(I,ceil(GraProx_spike_time/param.dt):end) = GraProx_spike_time;
         %%%%%%%%%%%%%%%%
         countrefracgraprox(I) = refracgraprox / param.dt; % neurons that just fired get
         % into the refractory period
         I = countrefracgraprox > 0;
         countrefracgraprox(I) = countrefracgraprox(I) - 1;
     
         I = find(countrefracgraprox == 0); % if countrefracgra = 0 the neuron can fire again
         % spike when V > thresh
         
         if ~isempty(I)
            aboveTh = Vgraprox(I,tt) >= Threshgraprox(I); %% [[Sam]] added spikelag
            countlag(aboveTh) = countlag(aboveTh) + 1;  %%[[Sam]] spike lag
            J = countlag == floor(spikelag/param.dt);
            if ~isempty(J)
              Vgraprox(I(J),tt) = param.SpikeV; % Action potential
              Sgraprox(I(J),tt) = 1; % Record spike time
              countlag(I(J)) = 0; %reset countlag for the spiking neurons;
            end
         end
         %%%% post V change effects, spike end %%%%
         
        
         %%%% recording intracellular currents %%%%
         if param.ProximalON == true 
         %now Iproxdist and Idistprox are not actively used in the calculation of membrane potential, but only for book keeping
         Iproxdist(:) = (Vgraprox(:,tt-1)-Vgradist(:,tt-1)) * GraDistal(1).wIntraDist; %wwIntraDist is the calculated, w is in the file; tt-1 instead of tt!!
         Iproxdist_matrix(:,tt-1) = Iproxdist(:); 
         
         Idistprox(:) = (Vgradist(:,tt-1)-Vgraprox(:,tt-1)) * GraProximal(1).wIntraProx; %wwIntraDist is the calculated, w is in the file
         Idistprox_matrix(:,tt)=Idistprox(:);
         
            
            
         end
         
         
         %%%% recording intracellular currents end %%%%
     end
     
    % Distal Granule cell potential
    % Forward Euler
%     Vgradist(:,tt) = Vgradist(:,tt - 1) + (param.dt ./ taugradist(:)) .* ...
%         ((ImitgradistAMPA(:) + ImitgradistNMDA(:) + (ImitgradistVDCC(:) - IVDCCBase) + Iproxdist(:))...
%         - Vgradist(:,tt - 1) + Restgradist(:) + Vnoisegradist);
    
    %[[Sam]] backwards Euler (partially, because only the compartmental coupling may cause stability problem, only that part is written out and solved for V at tt)
    
         if param.ProximalON == true
    Vgradist(:,tt) = tD*wD/(1+tD*wD+tP*wP) .* (Vgraprox(:,tt-1) + tP.*IP + ((1+tP*wP)/wD).*ID + ((1+tP*wP)/(tD*wD)*Vgradist(:,tt-1)));
         
         I = Vgraprox(:,tt-1)==param.SpikeV;
         Vgradist(I,tt) = Hypergraprox(I); %%%[[Sam]] setting Vgradist to Vhyper, can't set it earlier because then it would be changed again;

         else
             Vgradist(:,tt) = Vgradist(:,tt - 1) + (param.dt ./ taugradist(:)) .* ...
         ((ImitgradistAMPA(:) + ImitgradistNMDA(:) + (ImitgradistVDCC(:) - IVDCCBase) + Iproxdist(:))...
         - Vgradist(:,tt - 1) + Restgradist(:) + Vnoisegradist); 
         end
    % GABA release probability for Graded inhibition from distal gra dendrites
    Prelease_matrix(:,tt) = Prelease(CCa(:,tt),CCaBase,Threshgradist,CCaR);
    

end
    




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
% for debugging
%%%%

param.tAMPA0mit_gra = tAMPA0mit_gra;
param.t0pc_gra = t0pc_gra;
param.wMitGraProxAMPA = wMitGraProxAMPA;

% ORN and Mitral cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nMitral
    Mitral(ii).V = Vmit(ii,:); % Save neuronal activity
    Mitral(ii).S = Smit(ii,:); % Save spike time
end

% Granule cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nGraprox
     if param.ProximalON == true
         GraProximal(ii).V = Vgraprox(ii,:); % Save neuronal activity
         GraProximal(ii).S = Sgraprox(ii,:); % Save spike time
     else
         GraProximal(ii).V = zeros(1,round(param.tsim / param.dt)); % 0 vector
         GraProximal(ii).S = zeros(1,round(param.tsim / param.dt)); % 0 vector
     end
end
for ii = 1:param.nGradist
    if param.DistalON == true
        GraDistal(ii).V = Vgradist(ii,:); % Save neuronal activity
        GraDistal(ii).S = Sgradist(ii,:); % Save spike time
    else
        GraDistal(ii).V = zeros(1,round(param.tsim / param.dt)); % 0 vector
        GraDistal(ii).S = zeros(1,round(param.tsim / param.dt)); % 0 vector
    end
end


% Save input currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to proximal granule [[Sam]]
if param.ProximalON == true
    if param.DistalON == true
        InputCurrent.Idistprox = Idistprox_matrix;
    end
    
if param.CollateralON == true
    InputCurrent.ImitgraproxAMPA = ImitgraproxAMPA_matrix;
    InputCurrent.ImitgraproxNMDA = ImitgraproxNMDA_matrix; % added nmda
    InputCurrent.Vnoisegraprox = Vnoisegraprox_matrix;
    InputCurrent.Mg_block_graprox = Mg_block_graprox_matrix; %[Sam]
end
% pc to proximal granule [[Sam]]
if param.PCinputON == true
    InputCurrent.Ipcgraprox = Ipcgraprox_matrix;
    if strcmp(param.PCtype,'poisson')
        InputCurrent.PCspike = pcspikes_matrix;
    end
end

%    InputCurrent.Igragra = Igragra_matrix;
    InputCurrent.Igraspmit = Igraspmit_matrix;
end

% to distal granule
if param.DistalON == true
    if param.ProximalON == true
        InputCurrent.Iproxdist = Iproxdist_matrix;
    end
    InputCurrent.ImitgradistAMPA = ImitgradistAMPA_matrix;
    InputCurrent.ImitgradistNMDA = ImitgradistNMDA_matrix;
    InputCurrent.ImitgradistVDCC = ImitgradistVDCC_matrix;
    InputCurrent.Vnoisegradist = Vnoisegradist_matrix; %[Sam]
    InputCurrent.Imodgradist = Imodgradist_matrix; %[Sam]
    
    InputCurrent.Mg_block_gradist = Mg_block_gradist_matrix; %[Sam]
    
end


% to mitral
if param.DistalON == true
    InputCurrent.Igradistmit = Igradistmit_matrix;
%else
%    InputCurrent.Igraproxmit = Igraproxmit_matrix;
end
InputCurrent.Vnoisemit = Vnoisemit_matrix; %[[Sam]]
InputCurrent.Iext = Iext_matrix;

if param.PCinputON == true
    InputCurrent.Ipcgraprox = Ipcgraprox_matrix;
end

% Calcium currents
InputCurrent.Prelease = Prelease_matrix;
% InputCurrent.Pfiregraprox = Pfiregraprox;
InputCurrent.CCa = CCa;
InputCurrent.CCaBase = CCaBase;
InputCurrent.IVDCCBase = IVDCCBase;
InputCurrent.mgradist = m;
InputCurrent.hgradist = h;
InputCurrent.taum = taum_matrix; %[[Sam]]
toc
end



function Iext = SetExtInput(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function sets the external input to MCs
% The main function of this program is NeuroActivity.m
%
% Boleszek Osinski
% 03/05/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(2) % to be removed
Iext = zeros(param.nMitral,round(param.tsim / param.dt));
W = 0.001*(sort(rand(param.nMitral,1),'descend'))+param.Wmin; % hard-coded weights matrix
timerange = round(param.tinit / param.dt : round (param.tfinal / param.dt));
if param.flagRespiration == false %[[Sam]]
    Iext(:,timerange) = repmat(W .* (param.extparam + param.Inoise*randn(param.nMitral,1)),1,length(timerange)); %[[Sam]] added extparam; randn should be a vector!
else %[[Sam]]
    Iext(:,timerange) = (1/2 * param.extparam * sin(2*pi*param.RespFreq / (1/(param.dt/1000))* timerange) + param.extparam +param.Inoise*randn) .* W;
end
%rng shuffle %to be removed
end

function R = CreateRespFreq(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates an artificial respiration to modulate the bulb
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 04/06/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = param.dt:param.dt:param.tsim;
time = time / 1000; %converting the time to seconds
R = -cos(2 * pi * param.RespFreq * time);
R = R + 1;
R = R / 2;
end

function maxg = getmaxg(dt,tau1,tau2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the max amplitude for g, so we can use this value to
% normalize the curve
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/10/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 100;

x = 0:dt:tmax;
maxg = zeros(length(tau1),1);

for ii = 1: length(tau1)
    if ii == 1
        y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) * (exp(-(x)...
            / tau1(ii)) - exp(-(x) / tau2(ii)));
        maxg(ii) = max(y);
    else
        if tau1(ii) == tau1(ii -1) && tau2(ii) == tau2(ii -1)
            maxg(ii) = maxg(ii - 1);
        else
            y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii)))...
                * (exp(-(x) / tau1(ii)) - exp(-(x) / tau2(ii)));
            maxg(ii) = max(y);
        end
    end
end
end


function Ic = SetInoSpike_GraMit(gmax,P,E,V,W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns the distal granule summed graded input to mitral cells
% (weight matrix W is not necessarily square)
%
%
% Modifed by Boleszek Osinski on 07/16/2013 to allow for different numbers
% of neurons (particularly for nonspiking input from Gradist to Mitral)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: if W  = wGraMit then it has dimensions nMitral x nGranule

%     g = P .* gmax;
Ic = zeros(size(W,1),1);
    for ii = 1:size(W,1)
        
    Ic(ii) = sum(W(ii,:)' .* P * gmax(ii) * (E(ii) - V(ii))); % gmax is always 1 in practice
    end

    
end

function Ic = SetI(tau1,tau2,t,t0,normg,W,E,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the currents at each time step.
% This function is called from NeuroActivity.m
%
% Licurgo de Almeida
% 11/18/2010
%
% Ic: channel current
% gmax: maximum conductance
% tau1: channel's rising time
% tau2: channel's falling time
% t: step time
% t0: last time the pre synaptic neuron fired
% normg: normalizes the conductance curve between 0 and 1
% W: synaptic weights
% E: Nernst potential
% V: neuron's potential from last timestep
% Mcon: connection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g = (((tau1 * tau2) / (tau1 - tau2)) *...
    (exp(-(t - t0) / tau1) - exp(-(t - t0) / tau2))) / normg;

Ic = (W * g) .* (E - V);

end


function [INMDA,Mg_block] = SetI_NMDA(tau1,tau2,t,t0,gmax,W,E,V,param) %[Sam] added Mg_block output
% NMDA Synapse
g_norm = 1;
Mg_conc = param.Mg_conc;
gamma = param.gamma; % by default 0.016
eta = param.eta; %by default 0.28
E_Mg = 0; % assuming that Vrest = -70e-3

% % Code for visualizing all the variables of this current
% t1 = 2;
% t2 = 75;
% 
% t = 0:1:200;
% 
% v = (E_Mg-0.1):0.001:(E_Mg+0.1);
% 
% g = zeros(length(t), length(v));
% for i = 1:length(t)
%     for j = 1:length(v)
%         g(i,j) = g_norm * (exp(-t(i)/t2) - exp(-t(i)/t1))/(1+eta*Mg_conc*exp(-(v(j)-E_Mg)/gamma));
%     end
% end
% figure(4321)
% surf(v,t,g, 'FaceColor', 'interp', 'edgecolor', 'none', 'FaceLighting', 'phong')
% colorbar;
% camlight left;
% xlabel('V_M (V)');ylabel('t (ms)');zlabel('Conductance')
% 

% % variables for debugging
% tau1 = tauNMDA1(1);
% tau2 = tauNMDA2(1);
% t0 = tNMDA0mit_gra(:,tt);
% W = wMitGradist;
% E = ENMDAgradist(1);
% V = Vgradist(:,tt - 1);

% NOTE: If W is Wmitgradist then it has dimensions ngradist x nmit

% g = g_norm * (exp(-(t - t0)./tau1) - exp(-(t - t0)./tau2)); % dim nmit x 1
% Mg_block = 1./(1+eta*Mg_conc*exp(-(V-E_Mg)/gamma)); % dim ngradist x 1
% INMDA = (W * g) .* (V-E) .* Mg_block;

% Modified to allow for distributions over tau1 and tau2
Mg_block = 1./(1+eta*Mg_conc*exp(-(V-E_Mg)/gamma)); % dim ngradist x 1

INMDA = zeros(length(V),1);
for ii = 1:length(V)
    
    g = (exp(-(t - t0)./tau2(ii)) - exp(-(t - t0)./tau1(ii)))./gmax; % dim nmit x 1
    INMDA(ii) = sum(W(ii,:).* g') * (E-V(ii)) * Mg_block(ii);
end

% NOTE!!! All voltages are in V, not mV, so parameters such as gamma have
% to be scaled accordingly
end


function IVDCC = SetI_VDCC(W,E,V,m,h)
% L-type HVA Synapse
% defined as
% I_L = g*m*h*(V - E_L);

% % Visualize all variables of the current
% E_L = 0;
% V = 1e3*((E_L-0.1):0.001:(E_L+0.1)); % in mV
% 
% Compare L and N-type Ca models
% % Look at gating variable m
% % m' = (mbar - m)/taum;
% mbarL = 1./(1 + exp(-(V+50)./3));
% mbarN = 1./(1 + exp(-(V+45)./7));
% taumL = 18 * exp(-((V+45)./20).^2) + 1.5;
% taumN = 18 * exp(-((V+70)./25).^2) + 0.3;
% subplot(2,1,1)
% plot(V,mbarL,V,mbarN);xlim([-70 0])
% set(gca,'fontsize',14)
% legend('L-Type','N-Type','location','best')
% title('m')
% subplot(2,1,2)
% plot(V,taumL,V,taumN);xlim([-70 0])
% set(gca,'fontsize',14)
% title('\tau_m (ms)')
% xlabel('Vm_{GC} (mV)')
% 
% % plot N-type and NMDA together
% v = -0.1:0.001:0.1; % in V
% Mg_conc = 1;
% E_Mg = 0;
% gamma = 0.016;
% eta = 0.28;
% Mg_block_original = 1./(1+eta*Mg_conc*exp(-(v-E_Mg)/gamma));
% 
% scrsz = get(0,'ScreenSize');
% figH=figure;
% set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
% plot(V,mbarN,'-k',V,Mg_block_original,'.k');xlim([-80 0])
% set(gca,'fontsize',17)
% % legend('Bo''s model','NMDA')
% legend('N-Type','NMDA','location','southeast')
% legend boxoff
% xlabel('V_{rest,GC} (mV)');ylabel('Activation')


% 
% subplot(2,1,1)
% plot(V,mbar);xlim([-70 0])
% set(gca,'fontsize',17)
% xlabel('V (mV)');ylabel('Steady state m')
% subplot(2,1,2)
% plot(V,taum);xlim([-70 0])
% set(gca,'fontsize',17)
% xlabel('V (mV)');ylabel('\tau_{m} (ms)')
% 
% % Look at product of gating variables m*h
% Ccytmin = 0.1;
% Ccytmax = 0.5; % cytosolic [Ca] varries between 0.1 and 0.5 uM durring oscillations
% hmax = 0.00045/(0.00045+Ccytmax);
% hmin = 0.00045/(0.00045+Ccytmin);
% 
% plot(V,hmin*mbar,V,hmax*mbar);xlim([-70 0])
% set(gca,'fontsize',17)
% legend('m*h_{min}','m*h_{max}')
% legend boxoff
% xlabel('V (mV)');ylabel('Steady state m*h')

% Nernst potential for Ca
% RT = 300*8.31; % J/mole
% z= 2 ; % Ca ion valence
% F = 96485; % Faraday constant Coul/mole
% Cout = 1500;
% Cin = 0.5;
% 
% ECa = (RT/(z*F))*log(Cout/Cin)
% % ECa = 100 - 130 mV!


% I_L = g*m*h*(V - E_L)*s(t);
% NOTE: W has dimensions ngradist x nmit
% m has dimensions ngradist x 1

% Cin = 0.25; % cytosolic Ca concentration (uM)
% Cin = Cin(1);
% h = 0.00045/(0.00045+Cin);

IVDCC = zeros(length(V),1);
for ii = 1:length(V)
    IVDCC(ii) = sum(W(ii,:)) * m(ii) * h(ii) * (E(ii)-V(ii));
 %    IVDCC(ii) = max(W(ii,:)) * m(ii) * h(ii) * (E(ii)-V(ii));
end

% NOTE!!! All voltages are in V, not mV, so parameters have
% to be scaled accordingly
end

function P = Prelease(C,B,T,R)
% GABA release probablity
% 
% Boleslaw Osinski
% March 2015


% matrix of probabilities forced to be within range 0 - 1
% C - [Ca]
% B - baseline [Ca]
% T - threshold for maximum GABA release
P = (C - B - R) ./ (T - B - R);
% P = C ./ T;
J = P <= 0;
P(J) = 0;
J = P > 1;
P(J) = 1;

end




function Mat = SetConnections(cell1,cell2,cchance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the connection between different neurons
% The main function of this program is piriformmain.m
%
% Modified by Boleszek Osinski on 07/15/2013
%
% Licurgo de Almeida (original)
% 03/01/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: this function only sets fixed connections
        Mat = zeros(cell1,cell2);  %cell1/2 are integers indicating the number of cells
            for ii = 1:cell1
                convector = randperm(cell2);  
                convector = convector(1:round(cell2 * cchance)); %get a random subset of cells represented by cell2
                Mat(ii,convector) = 1; %set connection;
            end

end

function w = SetWeights(Mat,weights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the synaptic weights between different neurons
%
% % Modified by Boleszek Osinski on 07/15/2013
%
% Licurgo de Almeida (original)
% 03/02/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = zeros(size(Mat));
for ii = 1:length(weights)
    w(ii,:) = Mat(ii,:) * weights(ii);
end


end

%%%%%%%%%%%%%%%%Sam's functions%%%%%%%%%%%%%
%%%%set top down input to granule cells

function Ic = SetI_PC(param,tt)

Ic = zeros(param.nGraprox,round(param.tsim / param.dt));
W = 0.001*(sort(rand(param.nGraprox,1),'descend'))+param.Wmin;
graproxlist = 1:param.nGraprox;
timerange = round(param.tinit / param.dt : round (param.tfinal / param.dt));
phasedistribution = zeros(param.nGraprox,1);
if param.PCparam3 ~= 0
    unsync = datasample(graproxlist,param.nGraprox * param.PCparam3); %choose the ones receiving unsynchronous pcinput
    sync = graproxlist(find(ismember(graproxlist,unsync)));
    phasedistribution(unsync,:) = param.PCparam4 * randn(length(unsync),1);
end



if strcmp(param.PCtype,'sine')
    Ic(:,timerange) = (1/2 * param.PCparam1 * sin(2*pi*param.PCparam2 / (1/(param.dt/1000))* timerange + phasedistribution + 1.54) + param.PCparam1 +param.PCnoise*randn) .* W;
end


    
    
    

end
