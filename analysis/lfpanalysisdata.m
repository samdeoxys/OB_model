sampf = 2000;
data = odortrials;
location = 'OB1';
odor = 'emb';
time = 0; %odor onset at 1s
duration = 5; 
interval = [time*sampf+1:sampf*(time+duration)+1];
tobeanalyzed = data.(location).(odor)(:,interval)';
movingwin = [0.2 0.04];


params.tapers = [2,3];
params.fpass = [3,120];
params.Fs = sampf;
params.pad = 0;
params.trialave = 0;

% [S,f] = mtspectrumc(tobeanalyzed,params);
% figure
% hold on
% plot(f',S)
% title(['Power ',location,' ',odor,' ','time = ',num2str(time),'s ','duration = ',num2str(duration)])
% legend()
% set(gca,'YScale','log')
% hold off

[S,t,f] = mtspecgramc(tobeanalyzed(:,9),movingwin,params);
%[S,t,f] = mtspecgramc(tobeanalyzed(:,3),movingwin,params);
figure
pcolor(t,f,log(S(:,:)'))
title(['Power ',location,' ',odor,' ','time = ',num2str(time),'s ','duration = ',num2str(duration)])
colorbar
caxis([-8 -4])
colormap jet
shading interp
figure
plot(tobeanalyzed(:,9))