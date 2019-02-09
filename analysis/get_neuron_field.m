function mat = get_neuron_field(N,fieldname)
ncells = length(N);
npts = size(N(1).(fieldname),2);
mat = zeros(ncells,npts);

for i=1:ncells
    mat(i,:) = N(i).(fieldname);
end