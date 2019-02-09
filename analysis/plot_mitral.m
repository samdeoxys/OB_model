function h = plot_mitral(Mitral,InputCurrent,inds)
figure('un','norm','pos',[0.1,0.1,0.6,0.7])
hold on
L = [];
for i=1:length(inds)
plot(Mitral(inds(i)).V,'LineW',1.5)
plot(InputCurrent.Igradistmit(inds(i),:),'LineW',1.5)
plot(InputCurrent.Iext(inds(i),:),'LineW',1.5)
L = [L,strjoin(["Vmit",num2str(inds(i))]),strjoin(["Igaba",num2str(inds(i))]),strjoin(["Iext",num2str(inds(i))])];
end
legend(L)
datacursormode on

end