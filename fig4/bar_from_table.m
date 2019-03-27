function [b,ax] = bar_from_table(T,sliders,dataname,colors)
% grouped by frequency, left to right; b,g,r -> amplitude small to big;
% thickness of the color -> inputweight; 
data = T.(dataname);
%data = Tsmall.p_oscillation_driving_IPSP
data = reshape(data,[],length(sliders(3).value))'; % the order and ' are important!!! Reshaping by row, i.e. 100=4*25, with row 1 = 1-25, etc
%colors = {[0,0,1],[0,1,0],[1,0,0]};
%thickness = linspace(1,0.3,length((sliders(1).value)));
alphas = linspace(1,0.2,length(sliders(1).value));
figure
%hAx(1) = axes();
b = bar(data,'FaceColor','flat');
ax = gca;
set(gca,'XTickLabel',arrayfun(@(x)num2str(x),sliders(3).value,'uni',0))
title(dataname,'Interpreter','none')
ntrials = prod(arrayfun(@(i)length(sliders(i).value),1:length(sliders)));
for i =1:ntrials/length(sliders(3).value)
    ind_in_colors = find(sliders(1).value==T.(sliders(1).name)(i)); % map the amplitude of current trial to a color in colors
    ind_in_alpha = find(sliders(2).value ==T.(sliders(2).name)(i)); % map the inputweight of current trial to a thickness
    b(i).CData= repmat(colors{ind_in_colors},length(sliders(3).value),1);
    b(i).FaceAlpha = alphas(ind_in_alpha);
end
ind_sliders1legend =1:length(sliders(1).value); % indices for the graphs to display the legends for the parameter in sliders(1) 
sliders1legend = strsplit(num2str(sliders(1).value));
sliders1legend = cellfun(@(x)[sliders(1).name(end),' ',x],sliders1legend,'uni',0);
lgd = legend(b(ind_sliders1legend),sliders1legend);
% ind_sliders2legend = 2:length(sliders(1).value):1+length(sliders(1).value)*length(sliders(2).value);
% sliders2legend = strsplit(num2str(sliders(2).value));
% sliders2legend = cellfun(@(x)[sliders(2).name(end),' ',x],sliders2legend,'uni',0);
% lgd = legend(b([ind_sliders1legend,ind_sliders2legend]),[sliders1legend,sliders2legend]);
title(lgd,sliders(1).name(1:end-1))