function [figH,varargout] = plot_power(MitLFPs,GraDistLFPs,sampf,timevec)
% --Sam Zheng--
% Ploting the powerspectrum of LFPs simulated by IPSP on MC, EPSP on GCD,
% and V of MC
% ===========================
% Input:
% MitLFPs, GraDistLFPs = outputs of IandVLFP_GCE
% sampf, timevec = outputs of InitNetwork_GCE
% ===========================
% Output:
% figH = figure handle of the power spectrum
% varargout{1} = the FFT of IPSP on MC, optional output
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
gradFFTGI = fft(detrend(GraDistLFPs.MitGradistGlobal(trim:end-100),'constant'),NFFT)/L; 
mitFFTGV = fft(detrend(MitLFPs.VG(trim:end-100),'constant'),NFFT)/L;

% Plot %
scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
hold on
plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color','k','linestyle','-.'); 
%plot(f,2*abs(mitFFTEXTRA(1:NFFT/2+1)),'color','k','linestyle','--'); 
plot(f,2*abs(gradFFTGI(1:NFFT/2+1)),'color','k','linestyle',':'); 
plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color','k','linestyle','-');
plot(f,1e-4*ones(1,length(f)),'k--')

hold off
xlim(faxis);ylim([0 2e-3])
set(gca,'fontsize',fs)
legend('ILFP','ELFP','VmitLFP') %[Sam]
legend boxoff
xlabel('Frequency (Hz) ')        

varargout{1} = mitFFTGV;

end