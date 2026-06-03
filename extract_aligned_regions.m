% EXTRACT_ALIGNED_REGIONS  Parses a CLUSTAL alignment file to identify conserved residue positions.
%
% Input:
%   alignment_path - Path to a CLUSTAL .aln file produced by PyMOL
%
% Output:
%   aligned_regions - Row vector of residue indices (in the target protein) at positions
%                     marked conserved ('*') or similar ('.') in the alignment consensus line
%
% Reads the 'zero_prot' sequence row and its consensus line, removes gap positions
% from the target sequence, and returns the indices of structurally aligned residues.
function aligned_regions = extract_aligned_regions(alignment_path)
    % Read the CLUSTAL file
    fileID = fopen(alignment_path, 'r');
    tline = fgetl(fileID);
    line_counter = 1;
    elmnts_in_line = 63;
    first_char_in_line = 14;
    zero_prot_counter = nan;
    abs_start_index = 0;
    dash_counter = 0;
    aligned_regions = [];
    alignment_line = '';
    crnt_prot_line = '';
    while ischar(tline)
        if contains(tline, 'zero_prot')
            tline = fgetl(fileID);
            crnt_prot_line = [crnt_prot_line tline(first_char_in_line : end)];
            tline = fgetl(fileID);
            alignment_line = [alignment_line tline(first_char_in_line : end)];
        end
         tline = fgetl(fileID);      
    end
    crnt_prot_skips = regexp(crnt_prot_line, '-', 'all');
    alignment_line(crnt_prot_skips) = '';
    aligned_regions = regexp(alignment_line, '[*.]', 'all');
end
