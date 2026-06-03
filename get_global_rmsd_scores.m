% GET_GLOBAL_RMSD_SCORES  Computes whole-protein RMSD from CLUSTAL global alignment files.
%
% Input:
%   ds - Data structure with organism list and paths to global alignment files (_tx.aln)
%
% Output:
%   ds - Updated data structure; for each organism:
%        ds.crnt_<organism>_aligned_regions: conserved residue indices from the alignment
%        ds.crnt_<organism>_zero_global_rmsd_score: whole-protein RMSD in Angstrom
%
% RMSD is computed only over positions marked as conserved in the PyMOL CLUSTAL alignment.
function ds = get_global_rmsd_scores(ds)
    for j = 1: length(ds.organism_list)
        try
            crnt_organism = ds.organism_list{j};
            % Identify the residue positions over which RMSD is defined (from the global alignment).
            alignment_path = ds.(['crnt_' crnt_organism '_global_aln_path']);
            aligned_regions = extract_aligned_regions(alignment_path);
            ds.(['crnt_' crnt_organism '_aligned_regions']) = aligned_regions;
            %getting the RMSD
            crnt_dist_mat = ds.(['crnt_' crnt_organism '_zero_dist_matrix']);
            crnt_rmsd = get_rmsd_from_aa_array(crnt_dist_mat, aligned_regions);
            ds.(['crnt_' crnt_organism '_zero_global_rmsd_score']) = crnt_rmsd;
        catch Error
            continue
        end
    end
end