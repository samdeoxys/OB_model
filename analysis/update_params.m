function [Mitral,GraDistal,GraProximal,param]=update_params(Mitral,GraDistal,GraProximal,param,Tsliders)
%%% update parameters given the struct array; 
%%% naming rule for the fieldnames in Tsliders: eg 'Mitral_tau', the first
%%% part should be the name of the struct in which the parameter occurs,
%%% the second is the name of the parameter. 
%%% Tsliders  --- eg: Tsliders(1).Mitral_tau = 5, Tsliders(1).GraDistal_tau =
%%% 2;

nloops = height(Tsliders);
parnames = Tsliders.Properties.VariableNames;
npars = length(parnames);

for i =1:nloops
    for j=1:npars
        value = Tsliders.(parnames{j})(i);
        name = strsplit(parnames{j},'_'); % cell array of strings
        switch name{1}
            case 'Mitral'
                [Mitral.(name{2})] = deal(value);
            case 'GraDistal'
                [GraDistal.(name{2})] = deal(value);
            case 'GraProximal'
                [GraProximal.(name{2})] = deal(value);
            case 'param'
                param.(name{2}) = deal(value);
        end
    end
    
end
end

    
        
        
    