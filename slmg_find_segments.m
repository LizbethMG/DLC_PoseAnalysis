function segments = slmg_find_segments(nan_idx, gap_threshold)
        % Find continuous segments of non-NaN data with gaps smaller than gap_threshold
        segments = [];
        start_idx = 1;
        while start_idx <= length(nan_idx)
            if nan_idx(start_idx)
                start_idx = start_idx + 1;
                continue;
            end
            end_idx = start_idx;
            while end_idx < length(nan_idx) && sum(nan_idx(end_idx+1:min(end_idx+gap_threshold, length(nan_idx)))) < gap_threshold
                end_idx = end_idx + 1;
            end
            segments = [segments; start_idx, end_idx];
            start_idx = end_idx + 1;
        end
    end