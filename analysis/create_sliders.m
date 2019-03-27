function Tsliders = create_sliders(slider_array)
%%% create the struct array of sliders, the fieldname is the name of the parameter to be updated and the value is the value
%%% slider_array : the array of structs with 'name' and 'value' as fields,
%%% value is a vector of values and name is the name of the parameter for
%%% sweeping
%%%
%%% sliders : a struct array of N (=the total number of trials) structs,
%%% each field is the the parameter for sweep
nsliders = length(slider_array); % the number of types of parameters for sweep
combinations = slider_array(1).value; % initialize the combinations
for i=2:nsliders
    combinations = combvec(combinations,slider_array(i).value);
end

Tsliders = array2table(combinations','VariableNames',{slider_array.name});

    
end