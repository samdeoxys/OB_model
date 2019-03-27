function varargout = plot_power(MitLFPs,GraDistLFPs,sampf,timevec,varargin)
% --Sam Zheng--
% Ploting the powerspectrum of LFPs simulated by IPSP on MC, EPSP on GCD,
% and V of MC
% ===========================
% Input:
% MitLFPs, GraDistLFPs = outputs of IandVLFP_GCE
% sampf, timevec = outputs of InitNetwork_GCE
% varargin{1} = showplot = return the plot and handle if 1
% varargin{2} = plottitle = the title of the plot
% ===========================
% Output:
% 
% varargout{1} = the FFT of {IPSP on MC, Vgradist, Vmit}, optional output
% varargout{2} = NFFT, the number of points for FFT
% varargout{3} = f, the vector of frequencies
% varargout{4} = lfpnames
% varargout{5} = figure handle of the power spectrum

fs = 16;% fontsize
faxis = [15 120];
% params for FFT %
trim = 1000; % trim beginning and end to avoid edge effects
L = length(timevec(trim:end-100));  % Length of simulation
NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
f = sampf/2*linspace(0,1,NFFT/2+1);

% FFT %
mitFFTGI = fft(detrend(MitLFPs.GradistMitGlobal(trim:end-100),'constant'),NFFT)/L;
%mitFFTEXTRA = fft(detrend(MitLFPs.extra(trim:end-100),'constant'),NFFT)/L; 
%gradFFTGI = fft(detrend(GraDistLFPs.MitGradistGlobal(trim:end-100),'constant'),NFFT)/L;
gradFFTGV = fft(detrend(GraDistLFPs.VG(trim:end-100),'constant'),NFFT)/L;
mitFFTGV = fft(detrend(MitLFPs.VG(trim:end-100),'constant'),NFFT)/L;
LFPs = {mitFFTGI,gradFFTGV,mitFFTGV};
LFPs = cellfun(@(x)2*abs(x(1:NFFT/2+1)),LFPs,'uni',0);
linestyles = {'-.',':','-'};

% Plot %

if numel(varargin) == 0||(numel(varargin)>0&&varargin{1}==true)
    scrsz = get(0,'ScreenSize');
    figH=figure;
    set(figH,'position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
    hold on
    for i=1:length(LFPs)
        plot(f,LFPs{i},'color','k','linestyle',linestyles{i})
    end
    plot(f,1e-4*ones(1,length(f)),'k--')
    varargout{5} = figH;
    if numel(varargin) == 2
        title(varargin{2})
    end
end

hold off
xlim(faxis);ylim([0 2e-3])
set(gca,'fontsize',fs)
legend('ILFP','VgradistLFP','VmitLFP') %[Sam]
legend boxoff
xlabel('Frequency (Hz) ')        

varargout{1} = LFPs;
varargout{2} = NFFT;
varargout{3} = f;
varargout{4} = {'IPSP','Vgradist','Vmitral'};


end