% ADD_RATIO_COLUMN_TO_TABLES  Adds a global_rmsd_projections_by_global_rmsd ratio column
% to all coupling tables in a directory.
%
% One-time preprocessing step required before running evo_struct_align.m.
% The ratio (global_rmsd_projections / global_rmsd) is used by evo_struct_align.m
% to filter valid alignments (ratio < 2). The 28 tables in multi_evolution_tables/
% already have this column applied. Only re-run if tables are regenerated.
%
% Update folderPath below to point to multi_evolution_tables/ before running.
% Writes updated tables back to the same files in-place.
function add_ratio_column_to_tables(folderPath)
    % Get list of all .xlsx files in the folder
    % Update this path to the local multi_evolution_tables/ directory before running.
    folderPath = 'D:\multi_evolution\all_tables\';
    files = dir(fullfile(folderPath, '*.xlsx'));

    for k = 1:length(files)
        fileName = files(k).name;
        filePath = fullfile(folderPath, fileName);

        try
            % Read the table from Excel
            T = readtable(filePath);

            % Check if required columns exist
            if all(ismember({'global_rmsd_projections', 'global_rmsd'}, T.Properties.VariableNames))

                % Compute the new ratio column
                ratio = T.global_rmsd_projections ./ T.global_rmsd;

                % Replace 0/0 (which results in NaN) with 0
                ratio(isnan(ratio)) = 0;

                % Add the column to the table
                T.global_rmsd_projections_by_global_rmsd = ratio;

                % Write the updated table back to the same file
                writetable(T, filePath);

                fprintf('Processed: %s\n', fileName);
            else
                fprintf('Skipped (missing columns): %s\n', fileName);
            end

        catch ME
            fprintf('Error processing %s: %s\n', fileName, ME.message);
        end
    end
end