% find_covariance_peak: for two shifted (anti)correlated signals, s1 and s2,
% this fuction estimates the lag of s2 with respect to s1 by finding the
% maximum/minimum of the crosscovariance xcov(s1, s2)
function lag_samples = find_covariance_peak(covariance, max_lag, search_range, filter_length)

	% find maximum or minimum of filtered derivative
	dfdt = filtfilt(ones(100, 1)/100, 1, double(diff(covariance)));
	[dfdt_min dfdt_min_idx] = min(dfdt((max_lag + 1 - search_range):(max_lag + 1 + search_range)));
	[dfdt_max dfdt_max_idx] = max(dfdt((max_lag + 1 - search_range):(max_lag + 1 + search_range)));
	filtered_cov = filtfilt(hamming(filter_length)/sum(hamming(filter_length)), 1, double(covariance));
	if (dfdt_max_idx < dfdt_min_idx)
		% maximum of derivative is to the left -> covariance is positive
		[cov_max, idx] = max(filtered_cov((max_lag + 1 - search_range):(max_lag + 1 + search_range)));
	else
		% maximum of derivative is to the right -> covariance is negative
		[cov_min, idx] = min(filtered_cov((max_lag + 1 - search_range):(max_lag + 1 + search_range)));
	end
	lag_samples = -(idx - search_range - 1);

end
