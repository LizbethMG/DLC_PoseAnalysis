function vector_processed = slmg_interpolate_wo_gaps(vector, threshold_gap)

% slmg_interpolate_wo_gaps: Interpolates gaps in the input time series vectors x and y.
% The interpolation is performed only on gaps smaller than the specified threshold.
% Large gaps (>= threshold_gap) are left as NaN values.

% Inputs:
%   x - Time series vector with possible NaN values at the beginning, end, and within the data.
%   y - An additional time series vector with possible NaN values.
%   threshold_gap - Threshold for gap size; gaps smaller than this value will be interpolated.

% Outputs:
%   x_processed - The processed time series vector x with small gaps interpolated and large gaps left as NaN.
%   y_processed - The processed time series vector y with small gaps interpolated and large gaps left as NaN.


% PROCESS_VECTOR Helper function to process a single time series vector
% by interpolating gaps smaller than the specified threshold.

    % Remove leading and trailing NaNs
    start_idx = find(~isnan(vector), 1, 'first');
    end_idx = find(~isnan(vector), 1, 'last');
    vector_trimmed = vector(start_idx:end_idx);

    % Identify gaps
    % Identify gaps
    nan_idx = isnan(vector_trimmed);
    gap_starts = find(diff([0; nan_idx]) == 1);
    gap_ends = find(diff([nan_idx; 0]) == -1);
    gap_lengths = gap_ends - gap_starts + 1;

    % Create a copy for interpolation
    vector_interpolated = vector_trimmed;

 % Interpolate only gaps smaller than the threshold
    for i = 1:length(gap_lengths)
        if gap_lengths(i) < threshold_gap
            gap_start = gap_starts(i);
            gap_end = gap_ends(i);
            vector_interpolated(gap_start:gap_end) = interp1(find(~nan_idx), vector_trimmed(~nan_idx), gap_start:gap_end, 'spline');
        end
    end

    % Reconstruct the full vector with the original leading and trailing NaNs
    vector_processed = NaN(size(vector));
    vector_processed(start_idx:end_idx) = vector_interpolated;

end
