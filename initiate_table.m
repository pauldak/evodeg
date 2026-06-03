% INITIATE_TABLE  Creates an empty string table with column names from evolution_table_column_names.txt.
%
% Output:
%   evo_table - MATLAB table pre-allocated with 3500 rows of empty strings;
%               column names are read from evolution_table_column_names.txt
%               (one valid MATLAB identifier per line)
%
% Requires evolution_table_column_names.txt in the MATLAB working directory.
function evo_table = initiate_table()
    % Initialise table column names from external file
    evo_table_var_names = readlines(fullfile(fileparts(mfilename('fullpath')), 'data', 'evolution_table_column_names.txt'));
    evo_table_var_names = evo_table_var_names(strlength(evo_table_var_names) > 0);
    evo_table_var_names = matlab.lang.makeValidName(evo_table_var_names);
    numCols = numel(evo_table_var_names);
    numRows = 3500;  % Set to whatever number of rows you want
    data = repmat("", numRows, numCols);
    evo_table = array2table(data, 'VariableNames', evo_table_var_names);
end