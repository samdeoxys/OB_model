function [numhighext,varargout] = get_numhighext(Mitral,GraDistal,varargin)
% --Sam Zheng--
% For situations that involve having a proportion of GCs with high
% excitability and and the rest low excitability
% numhighext: nmit * 2 matrix: column 1 the number of high excitability GCs connected, column 2
% low excitability GCs connected
% varargin is Prelease, varargout is total prelease categorized into low
% ext group and high ext group
nmit = length(Mitral);
numhighext = zeros(nmit,2);
vrest_vec = get_neuron_field(GraDistal,'Vrest');
v_lowext = mode(vrest_vec);
for i=1:nmit
    inds = find(Mitral(i).ConnectionsGrad);
    numhighext(i,1) = sum(vrest_vec(inds)~= v_lowext);
    numhighext(i,2) = length(inds) - numhighext(i,1);
    if ~isempty(varargin)
        varargout{1}{i,1} = sum(varargin{1}(numhighext(i,1),:),1);
        varargout{1}{i,2} = sum(varargin{1}(numhighext(i,2),:),1);
    end
end

    
    
  