function c = heatmap_from_table(T,sliders,dataname,islog,samecolorsc)
% T -- the resulting data from triple_sweep
% sliders -- parameters used for sweeping
% dataname -- the dependent variable, as one variablename in T, to be
% ploted
% islog -- plot the log of the data if 1 by default, normal if 0
% samecolorsc -- same color scale across all plots if 1, default if 0

dims = arrayfun(@(x)length(sliders(x).value),1:3);
slicedtable = cell(dims(3),1);
%dataname = 'p_oscillation_driving_IPSP';
data = T.(dataname);
data = reshape(data,dims);
c=gobjects(dims(3),1);
for i =1:dims(3)
    figure('Name','heatmap')
    %pcolor(data(:,:,i))
    if islog==true
        imagesc(sliders(2).value,sliders(1).value,log(data(:,:,i)))
    else
        imagesc(sliders(2).value,sliders(1).value,data(:,:,i))
    end
    c(i)=gca;
    % fiddle with XDir, YDir reverse so that the ticklabels follow the
    % natural order; this depends on the order of numbers in the sliders
    set(gca,'YDir','normal')
    %shading interp
    title([sliders(1).name,' vs ',sliders(2).name,' ',sliders(3).name,'=',num2str(sliders(3).value(i))],'Interpreter','none')
    xlabel(sliders(2).name)
    ylabel(sliders(1).name)
    colorbar;
    %c(i).Label.String = 'power'
end
% below is for setting the colorbar to be the same across heatmaps; not
% helpful for within one heatmap
if samecolorsc ==true
    [cmin,cmax] = caxis(c(1));
    for i =1:dims(3)
        [ctempmin,ctempmax]=caxis(c(i));
        cmin = min(ctempmin,cmin);
        cmax = max(ctempmax,cmax);
    end
    for i=1:dims(3)
        caxis(c(i),[cmin,cmax])
    end
end
end