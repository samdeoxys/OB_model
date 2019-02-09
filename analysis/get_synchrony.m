function [corr_mean,varargout] = get_synchrony(N,tau,showfig)
trim = 1000;
Smat = get_neuron_field(N,'S');
convSmat = Smat;
corr_mat = zeros(length(N));
datarange = 1:10;
K = vRkernel(datarange,tau);
for i = 1:size(Smat,1)
    convSmat(i,:) = conv(Smat(i,:),K,'same');
end
corr_mat = corr(convSmat(:,trim:end)');
corr_mat = triu(corr_mat);
corr_vec = corr_mat(:);
corr_vec = corr_vec(corr_vec~=0 & corr_vec~=1);
corr_mean = mean(corr_vec(:),'omitnan');
if showfig 
    figure
    histogram(corr_vec(:),100)
    title('pairwise correlation of convolved spike train')

    hold on
    line([corr_mean,corr_mean],ylim,'LineWidth',2,'Color','r')
    hold off
end
varargout{1} = corr_mat;
end
    




function K = vRkernel(t,tau)
K = heaviside(t).*exp(-t./tau);

end