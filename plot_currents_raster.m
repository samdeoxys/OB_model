function plot_currents_raster(InputCurrent,Mitral,MitLFPs,ngradist,timevec,ntp,dt)
    cmat = ones(3)*1/ngradist; % bw
    fs = 16;% fontsize
    nmit = length(Mitral);
    vv=1; % not meaningful; 
    
    figure
    set(gcf,'position',[0,400,512,800]);
    xsc = [0 500];% x scale
    
    % plot currents
    ysc = [0 0.013]; % y scale
    % plot ImitgradistAMPA
    subplot(6,1,1)
    hold on
    for ii = 1:ngradist
%        plot(timevec,InputCurrent.ImitgradistAMPA(ii,:)/2,'color',ii*cmat(vv,:)); % color
         plot(timevec,InputCurrent.ImitgradistAMPA(ii,:)/2,'color',1-ii*cmat(vv,:)); % bw
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim(ysc)

    subplot(6,1,2)
    hold on
    for ii = 1:ngradist
       %plot(timevec,InputCurrent.ImitgradistNMDA(ii,:),'color',ii*cmat(vv,:)); % color
       plot(timevec,InputCurrent.ImitgradistNMDA(ii,:),'color',1-ii*cmat(vv,:)); % bw
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim(ysc)

    % plot ImitgradistVDCC
    subplot(6,1,3)
    hold on
    for ii = 1:ngradist
%        plot(timevec,InputCurrent.ImitgradistVDCC(ii,:),'color',ii*cmat(vv,:)); % color
        plot(timevec,InputCurrent.ImitgradistVDCC(ii,:),'color',1-ii*cmat(vv,:)); % bw
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim(ysc)

    % Prelease
    subplot(6,1,4)
    hold on
    for ii = 1:ngradist
%        plot(timevec,InputCurrent.Prelease(ii,:),'color',ii*cmat(vv,:)); % color
        plot(timevec,InputCurrent.Prelease(ii,:),'color',1-ii*cmat(vv,:))
    end
    hold off
    set(gca,'fontsize',fs)
    xlim(xsc);ylim([0 1])
    
    % Raster plot
    SpikeV = 65e-3;
        SPIKES = zeros(nmit,ntp);
        for n = 1:nmit
            SPIKES(n,:) = Mitral(n).S;
        end
    SPIKES = flipud(SPIKES); % order MCs in direction low E (bottom) to high E (top)
    subplot(6,1,5)
%     RasterPlot(SPIKES,dt,ntp*dt,ngradist*cmat(vv,:),fs,0); % color
    RasterPlot(SPIKES,dt,ntp*dt,'k',fs,0)
    xlim([xsc])

    % LFP plots
    subplot(6,1,6)
    hold on
%     plot(timevec,detrend(MitLFPs.GradistMitGlobal),'color',ngradist*cmat(vv,:)); % color
%     plot(timevec,detrend(MitLFPs.VG),'color',[ngradist*cmat(vv,:),0.2],'linewidth',2);
    plot(timevec,detrend(MitLFPs.GradistMitGlobal),'color','k','linewidth',2); % bw
    plot(timevec,detrend(MitLFPs.VG),'color',[0.6,0.6,0.6],'linewidth',2);
    hold off
    set(gca,'fontsize',fs)
    xlim([xsc]);ylim([-.005 0.005])
    xlabel('time (ms)')
     legend('ILFP','VLFP');legend boxoff
tightfig
end