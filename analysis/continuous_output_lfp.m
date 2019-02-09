a_list = [10 20 50];
c_list = -0.055;%[-0.06 -0.055 -0.05 -0.045 -0.04]
for j=1:3

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
                    
                    gradFFTGV = fft(detrend(GraDistLFPs.VG(trim:end-100),'constant'),NFFT)/L;
                    
                    sigmf = @(x,a,c)1./(1+exp(-a.*(x-c)));
                    
                    continuousout = zeros(size(GraDistal(1).V));
                    for i=1:720
                        continuousout = continuousout + sigmf(GraDistal(i).V,a_list(j),c_list(1));
                    end
                    continuousout = continuousout / 720;
                    continuousoutFFT = fft(detrend(continuousout(trim:end-100),'constant'),NFFT)/L;
                    
                    
                    
                    mitFFTGV = fft(detrend(MitLFPs.VG(trim:end-100),'constant'),NFFT)/L;
                    mitFFT1V = fft(detrend(MitLFPs.V1(trim:end-100),'constant'),NFFT)/L;
                    mitFFT2V = fft(detrend(MitLFPs.V2(trim:end-100),'constant'),NFFT)/L;
                    mitFFT3V = fft(detrend(MitLFPs.V3(trim:end-100),'constant'),NFFT)/L;

                    scrsz = get(0,'ScreenSize');
                    %figH1=figure;
                    %set(figH1,'position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
                    %hold off
                %     plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color',ngradist*cmat(vv,:)); % color
                %     plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color',[ngradist*cmat(vv,:),0.2],'linewidth',2)
                %     plot(f,2*abs(mitFFT1I(1:NFFT/2+1)),'color',ngradist*cmat(vv,:),'linestyle','--')
                %     plot(f,2*abs(mitFFT3I(1:NFFT/2+1)),'color',ngradist*cmat(vv,:),'linestyle',':','linewidth',2);
                    % plot(f,2*abs(mitFFT1V(1:NFFT/2+1)),f,2*abs(mitFFT3V(1:NFFT/2+1)));
                  %  plot(f,2*abs(mitFFTGI(1:NFFT/2+1)),'color','k'); % bw
                %   plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color',[0.75,0.75,0.75],'linewidth',2);
                   % plot(f,2*abs(mitFFT1I(1:NFFT/2+1)),'color','k','linestyle','-.');
                   % plot(f,2*abs(mitFFT3I(1:NFFT/2+1)),'color','k','linestyle',':','linewidth',2);

                    %plot(f,2*abs(mitFFTEXTRA(1:NFFT/2+1)),'color','k','linestyle','--'); %[[Sam]]
                   % plot(f,2*abs(gradFFTGI(1:NFFT/2+1)),'color','k','linestyle',':'); %[[Sam]]
                  %  plot(f,2*abs(mitFFTGV(1:NFFT/2+1)),'color','k','linestyle','-');
                   % plot(f,1e-4*ones(1,length(f)),'k--')
%                    figure('position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)])
%                    plot(f,2*abs(gradFFTGV(1:NFFT/2+1)),'color','k','linestyle','-')
%                     xlim(faxis(1,:));ylim([0 2e-3])
                    
                    
                   figure('position',[0,400,scrsz(3)-0.74*scrsz(3),scrsz(4)-0.7*scrsz(4)])
                   plot(f,2*abs(continuousoutFFT(1:NFFT/2+1)),'color','k','linestyle','-')
                   xlim(faxis(1,:));
                   title(num2str(a_list(j)))
end