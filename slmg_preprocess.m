function [x_smooth, y_smooth] = slmg_preprocess(x, y, zscore_threshold, window_size)

 % Remove Outliers Using Z-score Method
    x_mean = mean(x, 'omitnan');
    x_std = std(x, 'omitnan');
    y_mean = mean(y, 'omitnan');
    y_std = std(y, 'omitnan');

    x_zscore = (x - x_mean) ./ x_std;
    y_zscore = (y - y_mean) ./ y_std;

    x_outliers = abs(x_zscore) > zscore_threshold;
    y_outliers = abs(y_zscore) > zscore_threshold;

    x_clean = x;
    y_clean = y;
    x_clean(x_outliers) = NaN;
    y_clean(y_outliers) = NaN;

    % Plot the original and clean data for comparison
    figure;
    subplot(2,1,1);
    plot(x, 'r.'); hold on;
    plot(find(x_outliers), x(x_outliers), 'kx');
    legend('Original', 'Outliers');
    title('X Coordinates');
    
    subplot(2,1,2);
    plot(y, 'r.'); hold on;
    plot(find(y_outliers), y(y_outliers), 'kx');
    legend('Original', 'Outliers');
    title('Y Coordinates');

%     % Define the IQR multiplier threshold
%     iqr_multiplier = 1.5;
%     
%     % Calculate IQR for x and y coordinates
%     x_iqr = iqr(x);
%     y_iqr = iqr(y);
%     
%     % Calculate the quartiles
%     x_q1 = quantile(x, 0.25);
%     x_q3 = quantile(x, 0.75);
%     y_q1 = quantile(y, 0.25);
%     y_q3 = quantile(y, 0.75);
%     
%     % Define the bounds for outliers
%     x_lower_bound = x_q1 - iqr_multiplier * x_iqr;
%     x_upper_bound = x_q3 + iqr_multiplier * x_iqr;
%     y_lower_bound = y_q1 - iqr_multiplier * y_iqr;
%     y_upper_bound = y_q3 + iqr_multiplier * y_iqr;
%     
%     % Identify outliers
%     x_outliers = x < x_lower_bound | x > x_upper_bound;
%     y_outliers = y < y_lower_bound | y > y_upper_bound;
%     
%     % Replace outliers with NaN
%     x(x_outliers) = NaN;
%     y(y_outliers) = NaN;
%     % Define a gap threshold (number of consecutive NaNs to consider as a large gap)
%     gap_threshold = 10;
%     
%     % Identify segments with gaps smaller than the threshold
%     gap_indices = find(diff(find(isnan(x))) > gap_threshold);
%     segment_starts = [1; gap_indices + 1];
%     segment_ends = [gap_indices; length(x)];
    
%     % Initialize arrays for interpolated data
%     x_interp = x;
%     y_interp = y;
%     % Interpolate within segments
%     for i = 1:length(segment_starts)
%         start_idx = segment_starts(i);
%         end_idx = segment_ends(i);
%     
%         % Extract the segment
%         segment_x = x(start_idx:end_idx);
%         segment_y = y(start_idx:end_idx);
%         segment_frames = (start_idx:end_idx)';
%     
%         % Identify valid data points (non-NaN) within the segment
%         valid_indices_x = ~isnan(segment_x);
%         valid_indices_y = ~isnan(segment_y);
%     
%         % Apply spline interpolation within the segment if there are enough points
%         if sum(valid_indices_x) > 2
%             x_interp(start_idx:end_idx) = interp1(segment_frames(valid_indices_x), segment_x(valid_indices_x), segment_frames, 'spline');
%         end
%     
%         if sum(valid_indices_y) > 2
%             y_interp(start_idx:end_idx) = interp1(segment_frames(valid_indices_y), segment_y(valid_indices_y), segment_frames, 'spline');
%         end
%     end
%     
%     % Identify indices of large gaps and set them to NaN in the interpolated data
%     large_gap_indices = find(diff(find(isnan(x))) > gap_threshold);
%     for i = 1:length(large_gap_indices)
%         start_idx = find(isnan(x), large_gap_indices(i));
%         end_idx = find(isnan(x), large_gap_indices(i) + 1) - 1;
%         x_interp(start_idx:end_idx) = NaN;
%         y_interp(start_idx:end_idx) = NaN;
%     end
%     
%     % Apply a median filter to the interpolated data
%     window_size = 5; % Adjust the window size as needed
%     x_smooth = medfilt1(x_interp, window_size);
%     y_smooth = medfilt1(y_interp, window_size);
%     
%     % Plot the original, interpolated, and smoothed data for comparison
%     figure;
%     subplot(2,1,1);
%     plot(find(~isnan(x)), x(~isnan(x)), 'r.'); hold on;
%     plot(find(isnan(x)), x(isnan(x)), 'kx'); hold on;
%     plot(find(~isnan(x_interp)), x_interp(~isnan(x_interp)), 'b-');
%     plot(find(~isnan(x_smooth)), x_smooth(~isnan(x_smooth)), 'g-');
%     legend('Original', 'NaNs', 'Interpolated', 'Smoothed');
%     title('X Coordinates');
%     
%     subplot(2,1,2);
%     plot(find(~isnan(y)), y(~isnan(y)), 'r.'); hold on;
%     plot(find(isnan(y)), y(isnan(y)), 'kx'); hold on;
%     plot(find(~isnan(y_interp)), y_interp(~isnan(y_interp)), 'b-');
%     plot(find(~isnan(y_smooth)), y_smooth(~isnan(y_smooth)), 'g-');
%     legend('Original', 'NaNs', 'Interpolated', 'Smoothed');
%     title('Y Coordinates');
end