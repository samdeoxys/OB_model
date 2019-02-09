function varargout = analyze_firing(mit,fig_title)
trim=1000;
figure('un','norm','pos',[0.01,0.1,0.45,0.9])
Smit = get_neuron_field(mit,'S');
Smit_mean_vec = mean(Smit(:,trim:end),2)*10000;
Smit_std_vec = std(Smit(:,trim:end),0,2)*10000;
Smit_cv_vec = Smit_std_vec./Smit_mean_vec;
stat_Smit.mean = mean(Smit_mean_vec);
stat_Smit.std = mean(Smit_std_vec);
stat_Smit.cv = mean(Smit_cv_vec,'omitn');
Smit_master = {Smit_mean_vec,Smit_std_vec,Smit_cv_vec};
names = fieldnames(stat_Smit);
nbins = 15;
nfields = length(names);

for i=1:nfields
    subplot(nfields,1,i)
    histogram(Smit_master{i},nbins,'norm','prob')
    title(['Smit ',names{i},' ',fig_title])
    hold on
    line([stat_Smit.(names{i}),stat_Smit.(names{i})],get(gca,'YLim'),'color','r','LineS',':','LineW',1.0)
    legend('',num2str(stat_Smit.(names{i})))
    hold off
end
varargout{1} = stat_Smit;
tightfig
end