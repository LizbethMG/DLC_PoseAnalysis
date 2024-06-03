function [x_smooth, y_smooth, pre_proc_results] = slmg_preprocess(x, y, zscore_threshold, window_size, gap_threshold)
% slmg_preprocess - Remove outliers from DeepLabCut data using Z-score method,
% interpolate, and apply median filter while handling large gaps.
%
% Syntax: [x_smooth, y_smooth] = slmg_preprocess(x, y, zscore_threshold, window_size, gap_threshold)
%
% Inputs:
%    x - Original x coordinates.
%    y - Original y coordinates.
%    zscore_threshold - Z-score threshold for outlier detection.
%    window_size - Window size for the median filter.
%    gap_threshold - Threshold for defining large gaps in the data.
%
% Outputs:
%    x_smooth - Smoothed x coordinates.
%    y_smooth - Smoothed y coordinates.
%
% Example:
%    x = [data.bodypart_x];
%    y = [data.bodypart_y];
%    zscore_threshold = 3;
%    window_size = 5;
%    gap_threshold = 10;
%    [x_smooth, y_smooth] = slmg_preprocess(x, y, zscore_threshold, window_size, gap_threshold);

% Colors for plots:
light_blue = '#92DCE5';
raspberry = '#D81159';
quinacridone = '#8F2D56';
midnight_green = '#004E64';
xanthous = '#FFBC42';

%% 1. Remove outliers from DeepLabCut data using Z-score method.
fprintf (' >>> Removing outliers with %.2f  Z-score threshold for outlier detection... \n', zscore_threshold);

% Calculate the percentage of NaN values relative to the total number of
% elements from the original data
x_numNaNs = sum(isnan(x));% Count the number of NaN values
y_numNaNs = sum(isnan(y));
x_totalElements = length(x);% Calculate the total number of elements in the vector
y_totalElements = length(y);
x_percentageNaNs = (x_numNaNs / x_totalElements) * 100;% Calculate the percentage of NaN values
y_percentageNaNs = (y_numNaNs / y_totalElements) * 100;

% Calculate mean and standard deviation ignoring NaNs
x_mean = mean(x, 'omitnan');
x_std = std(x, 'omitnan');
y_mean = mean(y, 'omitnan');
y_std = std(y, 'omitnan');

% Calculate Z-scores
x_zscore = (x - x_mean) ./ x_std;
y_zscore = (y - y_mean) ./ y_std;

% Identify outliers based on Z-score threshold
x_outliers = abs(x_zscore) > zscore_threshold;
y_outliers = abs(y_zscore) > zscore_threshold;

% Replace outliers with NaN in both x and y (because it has to be
% consistent in the body part coordinates )
x_clean = x;
y_clean = y;
x_clean(x_outliers) = NaN;
x_clean(y_outliers) = NaN;
y_clean(y_outliers) = NaN;
y_clean(x_outliers) = NaN;

% Plot the original and cleaned data for comparison
figure;
subplot(2,1,1);
plot(x, 'Color', midnight_green); hold on;
plot(find(x_outliers), x(x_outliers), 'rx'); hold on;
plot(x_clean, 'Color', light_blue); hold on;
title('Original X Coordinates with Outliers Marked');
legend('Original',  'Outliers', 'Clean data');

subplot(2,1,2);
plot(y, 'Color', midnight_green); hold on;
plot(find(y_outliers), y(y_outliers), 'rx'); hold on;
plot(y_clean, 'Color', light_blue); hold on;
title('Original Y Coordinates with Outliers Marked');
legend('Original',  'Outliers', 'Clean data');

sgtitle('1.  Remove Outliers from data');

fprintf('     Numbers of outliers initially found : X coordinates - %d outliers, Y coordinates - %d outliers \n', sum(x_outliers ==1), sum(y_outliers == 1))
fprintf('     Numbers of NaN after removing outliers : X coordinates - %d outliers, Y coordinates - %d outliers \n', sum(isnan(x_clean)), sum(isnan(y_clean)))

%% 2. Interpolation: Apply spline interpolation while handling large gaps.

fprintf ('2. Apply spline interpolation while handling large gaps.\n')
fprintf('     Make evident NAN values.\n')
% First plot the cleaned data and mark nan values
[xnum_nans, ynum_nans] = slmg_markNAN(x_clean, y_clean);

fprintf('     Start appliying interpolation to fill small gaps.\n')
threshold_gap = 50;

x_interp = slmg_interpolate_wo_gaps(x_clean, threshold_gap);
y_interp = slmg_interpolate_wo_gaps(y_clean, threshold_gap);

% Plot the original and processed time series
figure;
subplot(2,1,1);
plot(x_interp, 'Color', raspberry); hold on;
plot(x_clean, 'Color', light_blue); hold on;
title('X Coordinates interpolated');
legend('Clean data', 'Interpolated');

subplot(2,1,2);
plot(y_interp, 'Color', raspberry); hold on;
plot(y_clean, 'Color', light_blue); hold on;
title('Y Coordinates interpolated');
legend('Clean data', 'Interpolated');

title('Interpolated Time Series');

%% 3. Median Filter: Smooth the interpolated data using a median filter.
window_size = 25; %nb_samples
x_smooth = medfilt1(x_interp, window_size);
y_smooth = medfilt1(y_interp, window_size);

% Plot the smoothed data
figure;
subplot(2,1,1);
plot(x_interp, 'Color', raspberry); hold on;
plot(x_smooth, 'Color', xanthous); hold on;
title('X Coordinates smoothed');
legend('Interpolated', 'Smoothed');

subplot(2,1,2);
plot(y_interp, 'Color', raspberry); hold on;
plot(y_smooth, 'Color', xanthous); hold on;
title('Y Coordinates smoothed');
legend('Interpolated', 'Smoothed');

title('SmoothedTime Series');

% Calculate the percentage of NaN values relative to the total number of
% elements from the original data

xs_numNaNs = sum(isnan(x_smooth));% Count the number of NaN values
ys_numNaNs = sum(isnan(y_smooth));
xs_totalElements = length(x_smooth);% Calculate the total number of elements in the vector
ys_totalElements = length(y_smooth);
xs_percentageNaNs = (xs_numNaNs / xs_totalElements) * 100;% Calculate the percentage of NaN values
ys_percentageNaNs = (ys_numNaNs / ys_totalElements) * 100;

%% 4. Recap of the pre-processing
% A few che points: 
if length(x) ~= length(x_smooth)
    error('The lengths of x and pre-processed x do not match.');
end

if x_percentageNaNs ~= y_percentageNaNs
    error('The % occlusion in x and y do not match.');
elseif  xs_percentageNaNs ~= ys_percentageNaNs
    error('The % occlusion in x pre-processes and y pre-processes do not match.');
end

% Level of noise before and after processing the original data 

% Logical index to ignore NaN values
validIndices = ~isnan(x) & ~isnan(x_smooth);
% Calculate Mean Squared Error (MSE), ignoring NaN values
mse_value = mean((x(validIndices) - x_smooth(validIndices)).^2);
% Calculate the standard deviation of the noise, ignoring NaN values
noise_std_before = std(x(validIndices) - mean(x(validIndices)));
noise_std_after = std(x_smooth(validIndices) - mean(x_smooth(validIndices)));
% Calculate Signal-to-Noise Ratio (SNR) improvement
signal_power = mean(x(validIndices).^2);
snr_before = signal_power / noise_std_before^2;
snr_after = signal_power / noise_std_after^2;
snr_improvement = 10 * log10(snr_after / snr_before);

% Convert time to minutes
fps = 25;
num_samples = length(x);
time_ms = (0:num_samples-1) * fps; % Each sample is 25 ms apart
time_minutes = time_ms / 60000;

fprintf ('______________________________\n')
fprintf ('Summary of the pre-processing done:\n')
fprintf ('    Recording duration: %.2f minutes at %.0f frames/seconds\n', time_minutes(end), fps)
fprintf('     Original data: \n');
fprintf('           Number of NaN values: %d\n', x_numNaNs);
fprintf('           Percentage of NaN values: %.2f%%\n', x_percentageNaNs);
fprintf('     After pre-processesing: \n');
fprintf('           Number of NaN values: %d\n', xs_numNaNs);
fprintf('           Percentage of NaN values: %.2f%%\n', ys_percentageNaNs);
fprintf('           Mean Squared Error (MSE) before and after smoothing: %.4f\n', mse_value);
fprintf('           Standard deviation of noise before smoothing: %.4f\n', noise_std_before);
fprintf('           Standard deviation of noise after smoothing: %.4f\n', noise_std_after);
fprintf('           SNR improvement: %.4f dB\n', snr_improvement);

pre_proc_results = [time_minutes(end), fps, x_percentageNaNs, ys_percentageNaNs, snr_improvement];
% Plot original vs pre-processed data 
figure;

subplot(2,1,1);
plot(time_minutes, x, 'Color', midnight_green); hold on;
plot(time_minutes, x_smooth, 'Color', xanthous); hold on;
title('X Coordinates smoothed');
legend('Interpolated', 'Smoothed');

subplot(2,1,2);
plot(time_minutes, y, 'Color', midnight_green); hold on;
plot(time_minutes, y_smooth, 'Color', xanthous); hold on;
title('Y Coordinates smoothed');
legend('Interpolated', 'Smoothed');

title('SmoothedTime Series');
end
