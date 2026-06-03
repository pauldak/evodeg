% INITIATE_VARIABLES  Creates and returns the pipeline data structure with default parameters.
%
% Output:
%   ds - Struct with the following fields:
%        prediction_threshold   - QCDPred score cutoff for degron classification (0.85)
%        default_chain          - PDB chain identifier ('A')
%        zero_organism          - Reference species name ('Saccharomyces_cerevisiae')
%        organism_list          - Cell array of all 7 species names
%        number_of_organisms    - Number of comparison species (7)
%        gap_coupling           - Interval boundary expansion applied after coupling (8 aa)
%        minimal_angstrom_dist  - Minimum Ca distance threshold for coupling search (2 A)
%        wrong_index_indicator  - Sentinel value for a failed coupling match (-1)
%        nan_interval           - String sentinel for a failed coupling ('-1_-1')
%        nan_alignment_score    - Numeric sentinel for a failed alignment (99)
%        minimal_pdb_size       - Minimum valid PDB file size in bytes (10000)
%        segment_size           - Sliding window length in amino acids (32)
%        qcdpred_penalty_matrix - Custom QCDPred-derived substitution matrix for localalign
function ds = initiate_variables()
    load('folder_paths.mat'); 
    ds.tmp_data_table = data_table;
    ds.prediction_threshold = 0.85;
    ds.default_chain = 'A';
    ds.zero_organism = 'Saccharomyces_cerevisiae';
    ds.organism_list = {'Saccharomyces_cerevisiae', 'Amphimedon_queenslandica', 'Caenorhabditis_elegans',...
                        'Nematostella_vectensis', 'Drosophila_melanogaster', 'Danio_rerio', 'Homo_sapiens'};
    ds.number_of_organisms = length(ds.organism_list);
    ds.gap_coupling = 8; %gap that will allow more range when matching intervals
    ds.minimal_angstrom_dist = 2; %minimal allowed distances between two aa so they are concidered matched as a couple
    ds.wrong_index_indicator = -1; %index indicating no coupling was found, probably as a major structural differences.
    ds.nan_interval = '-1_-1';
    ds.nan_alignment_score = 99;
    ds.minimal_pdb_size = 10000; %bytes
    ds.alignment_type = 0; % Reserved: 1 = valid alignment, 2 = optimized segment, 3 = segment not found.
    ds.segment_size = 32; %based on average degron size which was 31.5
    % Load the QCDPred-derived substitution matrix. Amino acid order: ARNDCQEGHILKMFPSTWYV,
    % consistent with the IUPAC ordering expected by MATLAB's localalign function.
    ds.qcdpred_penalty_matrix = readtable([main_folder 'docs\penalty_matrix_for_matlab.xlsx']);
    ds.qcdpred_penalty_matrix = table2array(ds.qcdpred_penalty_matrix);
end