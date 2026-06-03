% Computes the RMSD between a reference structure and a set of aligned residues in a target structure.
function rmsd = get_rmsd_from_aa_array(dist_mat, aa_indexes_tx)
    dist_values = [];
    for i = 1:length(aa_indexes_tx)
        crnt_aa_index = aa_indexes_tx(i);
        crnt_val = min(dist_mat(:, crnt_aa_index));
        dist_values = [dist_values; crnt_val];
    end
    rmsd = sqrt((sum(dist_values.^2))/length(dist_values));