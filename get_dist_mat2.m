% GET_DIST_MAT2  Computes the pairwise Ca-Ca distance matrix between two PDB structures.
%
% Inputs:
%   pdb_zero   - PDB structure for the reference protein (output of pdbread)
%   pdb_crnt   - PDB structure for the comparison protein (output of pdbread)
%   chain_zero - Chain identifier to use from pdb_zero (e.g. 'A')
%   chain_crnt - Chain identifier to use from pdb_crnt (e.g. 'A')
%
% Output:
%   dist_matrix - n x m matrix of Euclidean distances (Angstrom) between all Ca atoms
%                 of pdb_zero (n rows) and pdb_crnt (m columns)
function dist_matrix = get_dist_mat2(pdb_zero, pdb_crnt, chain_zero, chain_crnt)
    crnt_zero_CA_coords = get_carbon_alpha_matrix(pdb_zero, chain_zero);
    crnt_transx_CA_coords = get_carbon_alpha_matrix(pdb_crnt, chain_crnt);
    dist_matrix = calculate_CA_distances(crnt_zero_CA_coords, crnt_transx_CA_coords);
end

function CA_coords = get_carbon_alpha_matrix(PDB, chainID)
   % Initialize return variable
    CA_coords = [];
    % Check if chainID is provided; if not, use the first chain available
    if nargin < 3
        chainID = PDB.Model.Atom(1).chainID;
    end
    % Iterate through the Atom structure to find the CA of the specified amino acid
    for i = 1:length(PDB.Model.Atom)
        atom = PDB.Model.Atom(i);
        if strcmpi(atom.chainID, chainID) && strcmp(atom.AtomName, 'CA')
            % Extract coordinates
            CA_coords = [CA_coords; atom.X, atom.Y, atom.Z];
        end
    end
end

function dist_matrix = calculate_CA_distances(XYZ1, XYZ2)
    % XYZ1 is an n x 3 matrix of carbon alpha positions for the first set
    % XYZ2 is an m x 3 matrix of carbon alpha positions for the second set
    % Preallocate the distance matrix
    n = size(XYZ1, 1);
    m = size(XYZ2, 1);
    dist_matrix = zeros(n, m);
    % Calculate distances
    for i = 1:n
        for j = 1:m
            % Extract the coordinates
            ca1 = XYZ1(i, :);
            ca2 = XYZ2(j, :);
            % Calculate the Euclidean distance
            dist_matrix(i, j) = sqrt(sum((ca1 - ca2) .^ 2));
        end
    end
end