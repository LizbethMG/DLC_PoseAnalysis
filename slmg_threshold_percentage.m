function [percentages] = slmg_threshold_percentage (signal, threshold)

    % Count the total number of elements in the signal
    totalElements = numel(signal);
    
    % Count the number of NaN values in the signal
    nanCount = sum(isnan(signal));
    
    % Calculate the percentage of the signal that has NaNs
    nanPercentage = (nanCount / totalElements) * 100;
    
    % Find the elements that cross the threshold
    aboveThreshold = signal > threshold;
    belowThreshold = signal <= threshold;
    
    % Calculate the percentage of the signal that crosses the threshold
    aboveThresholdPercentage = (sum(aboveThreshold) / totalElements) * 100;
    
    % Calculate the percentage of the signal that does not cross the threshold
    belowThresholdPercentage = (sum(belowThreshold) / totalElements) * 100;
    
    % Groups the results in a single vector
    percentages = [belowThresholdPercentage, aboveThresholdPercentage, nanPercentage];
end

