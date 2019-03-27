directory = './sweepfigures/triplesweep/';
Tfilename = 'Excitability_vs_entrainment_125trials';
load(Tfilename)

%% bar plots
% grouped by frequency, left to right; b,g,r -> amplitude small to big;
% thickness of the color -> inputweight; 
dataname = 'p_oscillation_driving_IPSP';
colors = {[0 0 1],[0 1 1],[0 1 0],[1 1 0],[1 0 0]};
[b,ax] = bar_from_table(T,sliders,dataname,colors);

%% heat maps
% slice the table by the three sliders
dataname = 'p_oscillation_driving_IPSP';
islog = 1;
samecolorsc = 1;
c = heatmap_from_table(T,sliders,dataname,islog,samecolorsc); % heatmap_from_table(T,sliders,dataname,islog,samecolorsc)
for i = 1:length(c)
    xlabel(c(i),'Amplitude')
    ylabel(c(i),'Mean Input')
    ctitle = c(i).Title.String;
    ctitle = ['Mean Input vs Amplitude, f input = ',ctitle(end-1:end)];
    if islog==1
        ctitle = [ctitle,' log'];
    end
    c(i).Title.String = ctitle;
end
        

