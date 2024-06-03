function path = slmg_getExperimentPath(Experiments_list, animal, compound, dose, timepoint)
    % Initialize the path as an empty string
    path = '';
    
    % Iterate over the rows of the cell array
    for i = 1:size(Experiments_list, 1)
        % Check if the current row matches the selection criteria
        if Experiments_list{i, 1} == animal && strcmp(Experiments_list{i, 2}, compound) && ...
                Experiments_list{i, 3} == dose && Experiments_list{i, 4} == timepoint
            % Return the path if a match is found
            path = Experiments_list{i, 5};
            return;
        end
    end
    
    % If no match is found, display a warning message
    if isempty(path)
        warning('No matching experiment found for the given selection criteria.');
    end
end