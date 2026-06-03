% ALIGNMENT_MATLAB  Development script comparing BLOSUM62 vs QCDPred custom alignment scoring.
%
% Aligns two hardcoded example sequences using nwalign with two different scoring matrices:
%   1. Standard BLOSUM62
%   2. Custom QCDPred-derived substitution matrix (penalty_matrix_for_matlab.xlsx)
% Used to verify and validate the custom matrix before integrating it into the pipeline.
%
% Not part of the main pipeline. Update the readmatrix path to match the local file location.

seq_a = 'LDSWASLILQWFED';
seq_b = 'PEEWGKLIYQWVSR';


% Perform global alignment using BLOSUM64
[score1, alignment1] = nwalign(seq_a, seq_b, 'ScoringMatrix', blosum62);

% Display result
disp('Alignment using BLOSUM62:');
disp(alignment1);

%% 

mat_a = readmatrix("G:\My_Drive\LabData\Manuscripts\2024 Shir\figs\fig33\penalty_matrix_for_matlab.xlsx");  % Replace this with your actual matrix
alphabet = 'ARNDCQEGHILKMFPSTWYV';  % Make sure this matches mat_a

% Perform global alignment using custom matrix
[score2, alignment2] = nwalign(seq_a, seq_b, 'ScoringMatrix', mat_a);

% Display result
disp('Alignment using custom matrix mat_a:');
 disp(alignment2);

