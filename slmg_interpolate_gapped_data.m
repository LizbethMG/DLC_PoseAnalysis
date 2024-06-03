function interp_data = slmg_interpolate_gapped_data(data, gap_threshold)
        % Interpolate segments of data where gaps are smaller than gap_threshold
        interp_data = data;  % Initialize with the original data
        nan_idx = isnan(data);  % Find indices of NaNs
        segments = slmg_find_segments(nan_idx, gap_threshold);  % Find segments to interpolate

        for i = 1:size(segments, 1)
            segment = segments(i, :);
            interp_data(segment(1):segment(2)) = interp1(find(~nan_idx(segment(1):segment(2))), ...
                data(~nan_idx(segment(1):segment(2))), segment(1):segment(2), 'spline');
        end
    end