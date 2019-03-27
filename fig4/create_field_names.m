function fieldnames = create_field_names(name_cell)
idx = cellfun(@(x)1:numel(x),name_cell,'uni',0);
ind_combinations = combvec(idx{:}); %  nnames * ncombs
[nnames,ncombs] = size(ind_combinations);
out = arrayfun(@(x)name_cell{x}(ind_combinations(x,:))',(1:nnames),'uni',0); % 1 * nnames
out = horzcat(out{:}); % ncombs * nnames
fieldnames = arrayfun(@(x)strjoin(out(x,:),'_'),1:ncombs,'uni',0);

