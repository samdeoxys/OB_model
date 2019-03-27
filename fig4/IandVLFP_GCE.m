function [Mitral GraProximal GraDistal param InputCurrent MitLFPs GraDistLFPs] = IandVLFP_GCE(Mitral,GraProximal,GraDistal,param)

% Created by Boleslaw Osinski and modified by Sam Zheng
% Simulate LFP activity using membrane voltage / synaptic current

% INPUTS
%
% Mitral, GraProximal, GraDistal, param -- initialized structs of neurons
% and parameters output by OB_network_GCE.m

% OUTPUTS
%
% Mitral, GraProximal, GraDistal -- updated structs of neurons, with
% the time series of membrane potential and spiking activity filled in
% param -- the same param as input
% InputCurrent -- a struct of the time series of synaptic currents and other variables used in
% synaptic currents
% MitLFPs, GraDistLFPs -- structs of the time series of simulated LFP




%%% Must sum ALL currents into each neuron type
    [Mitral GraProximal GraDistal param InputCurrent] = NeuroActivity(Mitral,GraProximal,GraDistal,param);
    % Mitral LFP
    nmit  = length(Mitral);
    % Filter input currrents into Mit (ignore respiration)
    Gradistmit_filtered = zeros(nmit,length(Mitral(1).V));
    
    Mitral_filtered = zeros(nmit,length(Mitral(1).V));
    for n = 1:nmit
        Gradistmit_filtered(n,:) = InputCurrent.Igradistmit(n,:); %+ InputCurrent.Igraspmit(n,:); %[[Sam]] updated to include spike dependent gaba, no longer in use
        %Gradistmit_filtered(n,:) = smoothdata(Gradistmit_filtered(n,:),'movmean',50);
        Gradistmit_filtered(n,:) = smooth(Gradistmit_filtered(n,:),50);
        
        Mitral_filtered(n,:) = Mitral(n).V;
        %Mitral_filtered(n,:) = smoothdata(Mitral_filtered(n,:),'movmean',50); 
        Mitral_filtered(n,:) = smooth(Mitral_filtered(n,:),50); 
    end
    
    % ILFP
    MitLFPs.GradistMitGlobal = sum(Gradistmit_filtered,1)/nmit;
    if mod(nmit/3,1) == 0
    MitLFPs.GradistMit1 = sum(Gradistmit_filtered(1:(nmit/3),:),1)/(nmit/3);
    MitLFPs.GradistMit2 = sum(Gradistmit_filtered((1+nmit/3):2*(nmit/3),:),1)/(nmit/3);
    MitLFPs.GradistMit3 = sum(Gradistmit_filtered((1+2*nmit/3):end,:),1)/(nmit/3);
    end
    
    % VLFP
    MitLFPs.VG = sum(Mitral_filtered,1)/nmit;
    if mod(nmit/3,1) == 0
    MitLFPs.V1 = sum(Mitral_filtered(1:(nmit/3),:),1)/(nmit/3);
    MitLFPs.V2 = sum(Mitral_filtered((1+nmit/3):2*(nmit/3),:),1)/(nmit/3);
    MitLFPs.V3 = sum(Mitral_filtered((1+2*nmit/3):end,:),1)/(nmit/3);
    end
    
    % Distal Granule LFP
    ngra  = length(GraDistal);
    Gradist_filtered = zeros(ngra,length(GraDistal(1).V));
    for n = 1:ngra
        MitGradist_filtered(n,:) = (InputCurrent.ImitgradistAMPA(n,:)+InputCurrent.ImitgradistNMDA(n,:)+InputCurrent.ImitgradistVDCC(n,:)) / 3; %%[[Sam]], the sum of all synaptic currents into gradist
        %MitGradist_filtered(n,:) = smoothdata(MitGradist_filtered(n,:),'movmean',50);
        MitGradist_filtered(n,:) = smooth(MitGradist_filtered(n,:),50); % using smooth for compatibility with earlier versions, but need curve fitting toolbox
        
        Gradist_filtered(n,:) = GraDistal(n).V;
        %Gradist_filtered(n,:) = smoothdata(Gradist_filtered(n,:));   
        Gradist_filtered(n,:) = smooth(Gradist_filtered(n,:),50);
    end
    GraDistLFPs.VG = sum(Gradist_filtered,1)/ngra; %
    GraDistLFPs.MitGradistGlobal = sum(MitGradist_filtered,1) / ngra;%
    MitLFPs.extra = abs(GraDistLFPs.MitGradistGlobal) + abs(MitLFPs.GradistMitGlobal); %[[Sam]] the sum of both IPSC and EPSC;  
    

end

    