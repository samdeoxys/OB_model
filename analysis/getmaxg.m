function maxg = getmaxg(dt,tau1,tau2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the max amplitude for g, so we can use this value to
% normalize the curve
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/10/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 100;

x = 0:dt:tmax;
maxg = zeros(length(tau1),1);

for ii = 1: length(tau1)
    if ii == 1
        y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) * (exp(-(x)...
            / tau1(ii)) - exp(-(x) / tau2(ii)));
        maxg(ii) = max(y);
    else
        if tau1(ii) == tau1(ii -1) && tau2(ii) == tau2(ii -1)
            maxg(ii) = maxg(ii - 1);
        else
            y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii)))...
                * (exp(-(x) / tau1(ii)) - exp(-(x) / tau2(ii)));
            maxg(ii) = max(y);
        end
    end
end
end