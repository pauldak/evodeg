% ITERATOR_FOR_DEGRON_ALIGNMENT_COMPARISON  Sweeps the overlap threshold parameter
% to find the value that maximises the difference between real and random degron overlap.
%
% Tests overlap thresholds from 0.10 to 1.00 in steps of 0.05, calls
% degron_alignment_comparison at each threshold, and collects the resulting
% real-vs-random difference scores for inspection.
%
% Requires: degron_alignment_comparison on the MATLAB path and all its dependencies
% Produces: overlaps_arr and optimezed_diff_arr in the workspace; results printed to console

overlap_vals = 0.1:0.05:1.0;
overlaps_arr = [];
optimezed_diff_arr = [];
for i = 1:length(overlap_vals)
    overlap = overlap_vals(i);    
    % Call your main function
    optimezed_diff = degron_alignment_comparison(overlap);
    overlaps_arr = [overlaps_arr, overlap];
    optimezed_diff_arr = [optimezed_diff_arr, optimezed_diff];
    fprintf('Testing overlap threshold = %.2f resulted in %.4f\n', overlap, optimezed_diff);
end
a = [overlaps_arr', optimezed_diff_arr'];