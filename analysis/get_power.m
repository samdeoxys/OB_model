function varargout=get_power(MitLFPs,GraDistLFPs,sampf,timevec)
% varargout{1} = the FFT of {IPSP on MC, Vgradist, Vmit}, optional output
% varargout{2} = NFFT, the number of points for FFT
% varargout{3} = f, the vector of frequencies
% varargout{4} = lfpnames
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
varargout{1} = {2*abs(mitFFTGI(1:NFFT/2+1)),2*abs(gradFFTGV(1:NFFT/2+1)),2*abs(mitFFTGV(1:NFFT/2+1))};
varargout{2} = NFFT;
varargout{3} = f;
varargout{4} = {'IPSP_LFP','Vgradist_LFP','Vmitral_LFP'};
end