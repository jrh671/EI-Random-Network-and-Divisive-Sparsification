function [spike_rate_matrix] = SlidingAverage_s(spike_matrix,UpdatePeriod)

ms = (1e-3);
cs = (1e-2);
s = 1;

% Parameters
bin_duration = UpdatePeriod; 
target_bin_duration = s;
num_bins = round(target_bin_duration / bin_duration); % Number of bins in the target duration

% Pre-allocate the output matrix
num_rows = size(spike_matrix, 1);
num_columns = ceil(size(spike_matrix, 2) / num_bins);
spike_rate_matrix = zeros(num_rows, num_columns);

% Compute the spike rates
for row = 1:num_rows
    for col = 1:num_columns
        start_idx = (col - 1) * num_bins + 1;
        end_idx = min(col * num_bins, size(spike_matrix, 2));
        spike_rate_matrix(row, col) = sum(spike_matrix(row, start_idx:end_idx)) ;
    end
end
end