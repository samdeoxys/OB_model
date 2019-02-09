function [mean_high,mean_low] = high_low_compare(GraDistal,InputCurrent,field)
% use GraDistal to get the high/low indices
% S.(field) is the thing to be compared
vrest_vec = [GraDistal.Vrest];
low = mode(vrest_vec);
low_inds = find(vrest_vec==low);
high_inds = find(vrest_vec~=low);

mean_high = mean(InputCurrent.(field)(high_inds,:),1);
mean_low = mean(InputCurrent.(field)(low_inds,:),1);
figure('un','norm','pos',[0.1,0.1,0.6,0.7])
hold on
title(['high low comapre',field])
plot(mean_high,'LineW',1.5)
plot(mean_low,'LineW',1.5)
legend('high','low')
hold off

end

