function [Mitral GraProximal GraDistal param] = OB_network_GCE(inputFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sam Zheng 2018
%
% This is a modification of code originally developed by Licurgo de
% Almeida 2013 and modified by Boleslaw Osinski 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
% inputFile - txt file with parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS:
% Mitral, GraProximal, GraDistal - structures containing initialized neurons
% param      -  model parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting program

tic;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Reading parameters from input file and creating neurons
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param.inputFile = inputFile;
    
    % Open the input file for reading
    fid1 = fopen(param.inputFile,'r');
    if fid1 == -1
        msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
        return;
    end
    
    str = fgetl(fid1);
    
    % Get network parameters.
    
    while ~strcmpi(str,'neurons') %read what's before neurons
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % network parameters
            otherwise
                param = SetNetworkParameters(param,str);
        end  
        str = fgetl(fid1); %%read line by line while the file is open
    end
    
    
    
    % Create cells
    [param,Mitral,GraProximal,GraDistal] = CreateCells(param);
    %each return is a cell of structs, with value (array) and label

    str = fgetl(fid1);
    
    % Get cell parameters.
    while ~strcmpi(str,'end')
        
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % neuronal parameters
            case 'mitral'
                celltype = 'mitral'; %celltype is a mediator to signal reading the specific parameters
            case 'graproximal'
                celltype = 'graproximal';
            case 'gradistal'
                celltype = 'gradistal';
            otherwise
                switch celltype
                    case 'mitral'
                        Mitral = SetNeuronParameters(Mitral,param.nMitral,str);
                    case 'graproximal'
                        GraProximal = SetNeuronParameters(GraProximal,param.nGraprox,str);
                    case 'gradistal'
                        GraDistal = SetNeuronParameters(GraDistal,param.nGradist,str);    
                end
                
        end
        
        str = fgetl(fid1);
    end

    
    fclose(fid1); 
    fname = inputFile(1:end - 3); 
    fname = strcat(fname,'mat');
    save(fname,'Mitral','GraProximal','GraDistal','param'); 
end

function param = SetNetworkParameters(param,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----------OB--------------OB----------------OB--------------------------
%
% This function sets the parameters for the different neurons
% Modified by Sam Zheng 2018-2019
% Modified by Boleslaw Osinski on 06/14/2013
%
% Licurgo de Almeida
% 12/20/2010
%
% Information not related with the parameters of different neurons.
% Path: path where we save input and output file
% dt: timestep (in ms)
% tsim: simulation time (in ms)
% tinit: stimulus begin (in ms)
% tfinal: stimulus end (in ms)
% nMitral: number of mitral cells
% nGradist: number of granule distal synapses
% nGraprox: number of granule cell soma
% GraGracon = if true, granule cells connect to each other
% DistalON = if true graded inhibitory distal Granule dendrites are present
% ProximalON = if true proximal Granule soma are present
% BulbH = Bulb height (in distance units)
% BulbW = Bulb width (in distance units)
% NoiseMit = Mitral cell noise std
% NoiseGraprox = Granule proximal dendrite noise std
% NoiseGradist = Granule distal dendrite noise std
% hCaflag = if true, h (VDCC innactiavtion) depends on [Ca]
% ExFrac = fraction of excited dGCs
% Respiration = if true, the input is modulated by an oscillation
% representing the respiration
% RespFreq = Respiratory frequency
% Inoise = noise fraction of ORN input
% Wmin = minimum of MC input weight
% SpikeV = Spike voltage
% CChanceGraMit = Chance of connection between Gra and Mit
% CChanceGraGra = Chance of connection between Gra
% [modified by Sam Zheng]
% rhoCaVDCC = proportionality constant between [Ca] and IVDCC
% rhoCaNMDA = proportionality constant between [Ca] and INMDA
% PCinputON = if true GC soma receives top down input from PC
% PCtype = specifies the type of the PC input, 'constant'/'sine'/'poisson'
% PCparam1-4 = parameters for PCinput
% extparam = strength for external input on Mitral Cells
% PCnoise = noise for PCinput
% Mg_conc, eta, gamma = parameters for INMDA
% wMod = the proportion of leaky current blocked by modulation
% gMod = the synaptic weight of the modulation current on GC (temp)
% EMod = the reversal potential for modulation current on GC (temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    
    % OB parameters
    case 'path'
        param.outputPath = ParValue; %path
    case 'dt'
        param.dt = str2double(ParValue); %step
    case 'tsim'
        param.tsim = str2double(ParValue); %simulation time
    case 'tinit'
        param.tinit = str2double(ParValue); %stimulus begin
    case 'tfinal'
        param.tfinal = str2double(ParValue); %stimulus end
    case 'nmitral'
        param.nMitral = str2double(ParValue); %number of Mitral cells
    case 'ngradist'
        param.nGradist = str2double(ParValue); %number of granule distal synapses
    case 'ngraprox'
        param.nGraprox = str2double(ParValue); %number of Granule cell soma
    case 'gragracon'
        param.GraGracon = str2double(ParValue); %Gra-Gra connections
    case 'distalon'
        param.DistalON = str2num(ParValue); %flag distal gra dendrites
    case 'proximalon'
        param.ProximalON = str2num(ParValue); %flag proximal gra dendrites
    case 'bulbh'
        param.BulbH = str2double(ParValue); %Bulb height
    case 'bulbw'
        param.BulbW = str2double(ParValue); %Bulb width
    case 'noisemit'
        param.noisemit = str2double(ParValue); %noise Mitral
    case 'noisegraprox'
        param.noisegraprox = str2double(ParValue); %noise granule prox
    case 'noisegradist'
        param.noisegradist = str2double(ParValue); %noise granule dist
    case 'preset'
        param.flagpreset = str2num(ParValue); %flag preset positions
    case 'hcaflag'
        param.hCaflag = str2num(ParValue); %flag [Ca] dependence of h   
    case 'rhoca'
        param.rhoCa = str2num(ParValue); %proportionality factor between ICa and [Ca]    
    case 'exfrac'
        param.ExFrac = str2num(ParValue); %fraction of excited dGCs     
    case 'respiration'
        param.flagRespiration = str2num(ParValue); %flag respiratory modulation
    case 'respfreq'
        param.RespFreq = str2double(ParValue); %respiratory frequency
    case 'inoise'
        param.Inoise = str2double(ParValue); %noise fraction of ORN input
    case 'wmin'
        param.Wmin = str2double(ParValue); %minimum of MC EXT input weight
    case 'spikev'
        param.SpikeV = str2double(ParValue); %Spike voltage
    case 'cchancegramit'
        param.CChanceGraMit = str2double(ParValue); %chance of connection between Granule and Mitral cells
    case 'cchancemitgraprox'
        param.CChanceMitGraProx = str2double(ParValue); % [[Sam]] chance of connection between mitral cells and graprox;
    
    case 'cchancegragra'
        param.CChanceGraGra = str2double(ParValue); %chance of connection between Granule cells
    
    case 'pcinputon'
        param.PCinputON = str2num(ParValue); %%[[Sam]] PC input to GCs or not; use str2num for logical!!!vary
        
    case 'pcparam1'
        param.PCparam1 = str2double(ParValue); %%[[Sam]] PC param1, for SetI_PC, now constant current value (before weight)
    case 'pcparam2'
        param.PCparam2 = str2double(ParValue); %%[[Sam]] PC param1, for SetI_PC, now frequency 
    case 'collateralon'
        param.CollateralON = str2num(ParValue);
    case 'extparam'
        param.extparam = str2double(ParValue); %%[[Sam]] param for ext input to mitral cells
    case 'pctype'
        param.PCtype = ParValue; %% [[Sam]] type of PCinput
    case 'pcnoise'
        param.PCnoise = str2double(ParValue); %%[[Sam]] pc input noise
    case 'pcparam3'
        param.PCparam3 = str2double(ParValue); %%[[Sam]] PC param3, for SetI_PC, now for the proportion of graprox that receive input with different phases
    case 'pcparam4'
        param.PCparam4 = str2double(ParValue); %%[[Sam]] PC param4, for SetI_PC, now for the variance of distributino of the phase
    case 'mg_conc'
        param.Mg_conc = str2double(ParValue); %% [[Sam]] tweak to create mg free environment
    case 'eta'
        param.eta = str2double(ParValue); %%[[Sam]] tweak NMDA so that mg=1, nmda current is 1/4 of ampa, mg=0, nmda epsp is 6/5 of ampa
    case 'gamma'
        param.gamma = str2double(ParValue); %% [[Sam]] tweak NMDA
    case 'rhocavdcc'
        param.rhoCaVDCC = str2double(ParValue);
    case 'rhocanmda'
        param.rhoCaNMDA = str2double(ParValue);
        
    case 'wmod' %[Sam] neuromodulation
        param.wMod = str2double(ParValue);
    
    case 'gmod'  %[Sam] modulation by adding an inward current
       param.gMod = str2double(ParValue);
    
    case 'emod' %[Sam] modulation by adding an inward current
       param.EMod = str2double(ParValue);
    
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end

end

function [param,Mitral,GraProximal,GraDistal] = CreateCells(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function initiallizes the structure for each neuron type
%
% Modified by Boleslaw Osinski on 06/16/2013 and 08/13/2014
%
% Licurgo de Almeida
% 12/21/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    Mitral(param.nMitral).label = 'Mitral';
    for ii = 1:param.nMitral
        Mitral(ii).input = []; %no input for now
        Mitral(ii).label = 'Mitral';
    end
    
    GraProximal(param.nGraprox).label = 'GraProximal';
    GraDistal(param.nGradist).label ='GraDistal';
    for ii = 1:param.nGraprox
        GraProximal(ii).input = [];
        GraProximal(ii).label = 'GraProximal';
    end
    for ii = 1:param.nGradist
        GraDistal(ii).input = [];
        GraDistal(ii).label = 'GraDistal';
    end


end

function N = SetNeuronParameters(N,ncells,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
% Modified by Sam Zheng 2018-2019
%
%
%
% * tau = charging time constant of the neuron (ms). ms (not s) is the basic time unit in this program.
% * Fthresh = Firing threshold in V
% * CCaTh = Threshold [Ca] for maximum GABA release
% * Vrest = resting potential
% * Vhyper = hyperpolarization potential
% * EAMPA = AMPA's Nernst potential.
% * ECa = Ca Nernst potential.
% * EGABA = GABA's Nernst potential.
% * tauAMPA1 = AMPA's rising tau.
% * tauAMPA2 = AMPA's falling tau.
% * tauGABA1 = GABA's rising tau.
% * tauGABA2 = GABA's falling tau.
%
%
% * gmaxAMPA = AMPA's max conductance
% * gmaxGABA = GABA's max conductance
% * gmaxGABAP = Periglomerula GABA's max conductance (for Mitral cells
% only)
% * IACh = Addition current when ACh is ON in non-spiking cells.
% * nGracon = number of granule cells connected to mitral or granule cells
% * tauPROX1 = proximal Mit-Gra synapse rising tau (if distal dendrites removed)
% * tauPROX2 = proximal Mit-Gra synapse falling tau (if distal dendrites removed)
% * tauAMPA1 = distal Mit-Gra AMPA receptor rising tau (Gra cell only)
% * tauAMPA2 = distal Mit-Gra AMPA receptor falling tau (Gra cell only)
% * tauNMDA1 = distal Mit-Gra NMDA receptor rising tau (Gra cell only)
% * tauNMDA2 = distal Mit-Gra NMDA receptor falling tau (Gra cell only)
% * tauVDCC = time constant of VDCC activation variable
% * wAMPAMI = excitatory synaptic weight from Mitral cell to Granule cell AMPA
% * wNMDAMI = excitatory synaptic weight from Mitral cell to Granule cell NMDA
% * wVDCCMI = synaptic conductance of Granule cell VDCC
% * wGABAGR = inhibitory synaptic weight from Granule cell to Gra or Mit
% * wIntraDist = synaptic weight from graproximal to gradistal [[Sam]]
% * wIntraProx = synaptic weight from gradistal to graproximal [[Sam]]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    % OB and PC
    case 'tau'
        for ii = 1:ncells
            N(ii).tau = str2double(ParValue); %time in ms
        end
    case 'fthresh'
        for ii = 1:ncells
            N(ii).FThresh = str2double(ParValue); %threshold in V
        end        
    case 'ccath'
        for ii = 1:ncells
            N(ii).CCaTh = str2double(ParValue); %[Ca] in uM
        end
    case 'vrest'
        for ii = 1:ncells
            N(ii).Vrest = str2double(ParValue); %potential in volts
        end
    case 'vhyper'
        for ii = 1:ncells
            N(ii).Vhyper = str2double(ParValue); %potential in volts
        end
    case 'noise'
        for ii = 1:ncells
            N(ii).Noise = str2double(ParValue);
        end
    case 'eampa'
        for ii = 1:ncells
            N(ii).EAMPA = str2double(ParValue); %AMPA (Na) reversal potential in volts
        end
    case 'eca'
        for ii = 1:ncells
            N(ii).ECa = str2double(ParValue); %Ca reversal potential in volts
        end
    case 'egaba'
        for ii = 1:ncells
            N(ii).EGABA = str2double(ParValue); %potential in volts
        end
    case 'tauampa1'
        for ii = 1:ncells
            N(ii).tauAMPA1 = str2double(ParValue); %time in ms
        end
    case 'tauampa2'
        for ii = 1:ncells
            N(ii).tauAMPA2 = str2double(ParValue); %time in ms
        end
    case 'taugaba1'
        for ii = 1:ncells
            N(ii).tauGABA1 = str2double(ParValue); %time in ms
        end
    case 'taugaba2'
        for ii = 1:ncells
            N(ii).tauGABA2 = str2double(ParValue); %time in ms
        end
        
        % OB
    case 'gmaxampa'
        for ii = 1:ncells
            N(ii).gmaxAMPA = str2double(ParValue); %AMPA channel
            % conductance in siemens
        end
    case 'gmaxgaba'
        for ii = 1:ncells
            N(ii).gmaxGABA = str2double(ParValue); %GABA channel
            % conductance in siemens
        end
    case 'ngracon'
        for ii = 1:ncells
            N(ii).NGraCon = str2double(ParValue); %number of granule cells
            % connected to a mitral or granule cell
        end
    case 'tauprox1'
        for ii = 1:ncells
            N(ii).tauPROX1 = str2double(ParValue); %time in ms
        end
    case 'tauprox2'
        for ii = 1:ncells
            N(ii).tauPROX2 = str2double(ParValue); %time in ms
        end
    case 'taudist1'
        for ii = 1:ncells
            N(ii).tauAMPA1 = str2double(ParValue); %time in ms
        end
    case 'taudist2'
        for ii = 1:ncells
            N(ii).tauAMPA2 = str2double(ParValue); %time in ms
        end
    case 'taunmda1'
        for ii = 1:ncells
            N(ii).tauNMDA1 = str2double(ParValue); %time in ms
        end
    case 'taunmda2'
        for ii = 1:ncells
            N(ii).tauNMDA2 = str2double(ParValue); %time in ms
        end
    case 'tauvdcc'
        for ii = 1:ncells
            N(ii).tauVDCC = str2double(ParValue); %time in ms
        end        
    case 'wampami'
        for ii = 1:ncells
            N(ii).wAMPAMI = str2double(ParValue); % excitatory synaptic weight from Mitral to Granule AMPA
        end
    case 'wnmdami'
        for ii = 1:ncells
            N(ii).wNMDAMI = str2double(ParValue); % excitatory synaptic weight from Mitral to Granule NMDA
        end
     case 'wvdccmi'
        for ii = 1:ncells
            N(ii).wVDCCMI = str2double(ParValue); % VDCC synaptic conductance
        end
    case 'wgabagr'
        for ii = 1:ncells
            N(ii).wGABAGR = str2double(ParValue); % inhibitory synaptic weight from Granule to Gra or Mit
            
        end
        % New parameters must be added here...
    case 'wintradist'
        for ii = 1:ncells
            N(ii).wIntraDist = str2double(ParValue); % [[Sam]]
        end
    case 'wintraprox'
        for ii = 1:ncells
            N(ii).wIntraProx = str2double(ParValue); %[[Sam]]
        end
        

        
    case 'tauca'
        for ii=1:ncells
            N(ii).tauCa = str2double(ParValue);
        end
    %[[Sam]] Not uesful in the current version
    case 'wgabagrsp'
        for ii = 1:ncells
            N(ii).wGABAGRSP = str2double(ParValue); 
        end
        
    case 'ccar'
        for ii=1:ncells
            N(ii).CCaR = str2double(ParValue);
        end
    case 'spikelag'
        for ii=1:ncells
            N(ii).spikelag = str2double(ParValue);
        end
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end

end
