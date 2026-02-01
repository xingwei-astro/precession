function [pks, locs, w, p] = findpeaks(data, varargin)
%FINDPEAKS Find local peaks in data
%   [PEAKS, LOCS] = FINDPEAKS(DATA) returns local maxima (peaks) in DATA
%   and their locations (indices)
%
%   [PEAKS, LOCS, WIDTHS, PROMINENCES] = FINDPEAKS(DATA) also returns
%   peak widths and prominences
%
%   FINDPEAKS(..., 'MINPEAKHEIGHT', MPH) finds peaks greater than MPH
%   FINDPEAKS(..., 'MINPEAKDISTANCE', MPD) sets minimum distance between peaks
%   FINDPEAKS(..., 'THRESHOLD', TH) sets minimum difference from neighbors
%   FINDPEAKS(..., 'NPEAKS', N) returns only N highest peaks
%   FINDPEAKS(..., 'SORTSTR', 'ascend'|'descend') sorts peaks
%
%   Example:
%       x = linspace(0, 10, 100);
%       y = exp(0.2*x) .* sin(3*x) + 0.1*randn(1,100);
%       [pks, locs] = findpeaks(y, 'MinPeakDistance', 20);
%       plot(x, y, 'b-', x(locs), pks, 'ro');

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'data', @isnumeric);
    addParameter(p, 'MinPeakHeight', -Inf, @isnumeric);
    addParameter(p, 'MinPeakDistance', 1, @isnumeric);
    addParameter(p, 'Threshold', 0, @isnumeric);
    addParameter(p, 'NPeaks', Inf, @isnumeric);
    addParameter(p, 'SortStr', 'descend', @ischar);
    parse(p, data, varargin{:});
    
    params = p.Results;
    
    % Initialize outputs
    pks = [];
    locs = [];
    w = [];
    p = [];
    
    % Ensure data is a column vector
    data = data(:);
    n = length(data);
    
    % Find all local maxima (simple approach)
    all_locs = find_local_maxima(data, params.Threshold);
    
    if isempty(all_locs)
        return;
    end
    
    % Get peak values
    all_pks = data(all_locs);
    
    % Apply MinPeakHeight filter
    height_mask = all_pks >= params.MinPeakHeight;
    all_locs = all_locs(height_mask);
    all_pks = all_pks(height_mask);
    
    if isempty(all_locs)
        return;
    end
    
    % Apply MinPeakDistance filter
    [all_locs, all_pks] = apply_min_peak_distance(all_locs, all_pks, params.MinPeakDistance);
    
    % Sort peaks
    [all_pks, sort_idx] = sort(all_pks, params.SortStr);
    all_locs = all_locs(sort_idx);
    
    % Apply NPeaks filter
    n_peaks = min(params.NPeaks, length(all_pks));
    pks = all_pks(1:n_peaks);
    locs = all_locs(1:n_peaks);
    
    % Calculate peak prominences if requested
    if nargout >= 4
        p = calculate_prominences(data, locs);
    end
    
    % Calculate peak widths if requested
    if nargout >= 3
        w = calculate_peak_widths(data, locs);
    end
    
    % Sort by location if not already
    [locs, sort_idx] = sort(locs);
    pks = pks(sort_idx);
    if ~isempty(w), w = w(sort_idx); end
    if ~isempty(p), p = p(sort_idx); end
end

%% Helper Functions

function locs = find_local_maxima(data, threshold)
    % Find indices where data is greater than both neighbors
    n = length(data);
    locs = [];
    
    for i = 2:n-1
        if data(i) > data(i-1) && data(i) > data(i+1) && ...
           (data(i) - data(i-1)) >= threshold && ...
           (data(i) - data(i+1)) >= threshold
            locs = [locs; i];
        end
    end
end

function [locs, pks] = apply_min_peak_distance(locs, pks, min_dist)
    % Sort by peak height (descending)
    [pks_sorted, sort_idx] = sort(pks, 'descend');
    locs_sorted = locs(sort_idx);
    
    % Keep peaks that are at least min_dist apart
    keep = true(size(locs_sorted));
    
    for i = 1:length(locs_sorted)
        if keep(i)
            % Find peaks too close to this one
            too_close = abs(locs_sorted - locs_sorted(i)) < min_dist;
            too_close(i) = false;  % Don't mark current peak
            
            % Mark them for removal (except the current one)
            keep(too_close) = false;
        end
    end
    
    locs = locs_sorted(keep);
    pks = pks_sorted(keep);
end

function prominences = calculate_prominences(data, locs)
    % Calculate peak prominence (simplified version)
    n = length(locs);
    prominences = zeros(n, 1);
    
    for i = 1:n
        peak_loc = locs(i);
        peak_val = data(peak_loc);
        
        % Find left and right bases
        left_base = find_left_base(data, peak_loc);
        right_base = find_right_base(data, peak_loc);
        
        % Prominence is height above the higher base
        base_height = max(data(left_base), data(right_base));
        prominences(i) = peak_val - base_height;
    end
end

function left_base = find_left_base(data, peak_loc)
    % Find lowest point to the left before a higher point
    peak_val = data(peak_loc);
    left_base = peak_loc;
    
    for i = peak_loc-1:-1:1
        if data(i) > peak_val
            break;
        end
        if data(i) < data(left_base)
            left_base = i;
        end
    end
end

function right_base = find_right_base(data, peak_loc)
    % Find lowest point to the right before a higher point
    peak_val = data(peak_loc);
    right_base = peak_loc;
    
    for i = peak_loc+1:length(data)
        if data(i) > peak_val
            break;
        end
        if data(i) < data(right_base)
            right_base = i;
        end
    end
end

function widths = calculate_peak_widths(data, locs)
    % Calculate peak width at half prominence
    n = length(locs);
    widths = zeros(n, 1);
    
    for i = 1:n
        peak_loc = locs(i);
        peak_val = data(peak_loc);
        
        % Calculate half-height
        left_base = find_left_base(data, peak_loc);
        right_base = find_right_base(data, peak_loc);
        base_height = max(data(left_base), data(right_base));
        half_height = base_height + (peak_val - base_height) / 2;
        
        % Find left and right intersections with half-height
        left_width = find_width_intercept(data, peak_loc, half_height, 'left');
        right_width = find_width_intercept(data, peak_loc, half_height, 'right');
        
        widths(i) = right_width - left_width;
    end
end

function intercept = find_width_intercept(data, peak_loc, half_height, direction)
    % Find where data crosses half_height
    if strcmp(direction, 'left')
        indices = peak_loc-1:-1:1;
    else
        indices = peak_loc+1:length(data);
    end
    
    intercept = peak_loc;
    
    for i = indices
        if (direction == 'left' && data(i) <= half_height) || ...
           (direction == 'right' && data(i) <= half_height)
            % Linear interpolation for more accurate intercept
            if i > 1 && i < length(data)
                if direction == 'left'
                    x1 = i; y1 = data(i);
                    x2 = i+1; y2 = data(i+1);
                else
                    x1 = i-1; y1 = data(i-1);
                    x2 = i; y2 = data(i);
                end
                intercept = x1 + (half_height - y1) * (x2 - x1) / (y2 - y1);
            else
                intercept = i;
            end
            break;
        end
    end
end