% load parameters
innFLUX_parameters;

results = load([output_folder '/results.mat']);
for (i = 1:length(results.TRACER))
	results.TRACER(i).lagtime2 = NaN(length(results.time), 1);
	results.TRACER(i).flux2 = NaN(length(results.time), 1);
end

% create folder for figures of cumulated covariances
if (exist([output_folder '/figures'], 'dir') == 0)
    if (mkdir([output_folder '/figures']) == 0)
        fprintf('could not create folder for figures: %s\n', [output_folder '/figures']);
        return;
    end
end

% cumulated covariances for data segments
if (isempty(data_segments))
    data_segments = [results.time(1) results.time(end)];
end
cumul_cov = zeros(length(results.TRACER), 2*MAX_LAG+1, length(data_segments)-1);
cumul_abs_cov = zeros(length(results.TRACER), 2*MAX_LAG+1, length(data_segments)-1);

intervals_per_day = int32(floor(86400*SAMPLING_RATE_SONIC/WINDOW_LENGTH));

% list covariance files
cov_files = list_files_datestring([output_folder '/cov'], {'cov'}, '.mat', 'yyyy.mm.dd');
fprintf('innFLUX_step2: found %d covariance files\n', length(cov_files));

% loop over all covariance files and cumulate covariances
fprintf('step 1 of 3: cumulating covariances\n');
for (n = 1:length(cov_files))

	% load covariance file
	cov_data = load([cov_files(n).path '/' cov_files(n).name]);
	
	% loop over all tracers
	for (i = 1:length(cov_data.TRACER))
		
        if (~isempty(tracer_indices) && sum(i == tracer_indices) == 0)
            % skip this tracer
            continue;
        end
		
		% find position in results file corresponding to current covariance file
		[~, results_idx] = min(abs(cov_files(n).date - results.time));
		
		% cumulate covariances
		for (k = 1:intervals_per_day)
            linenum = results_idx-1+k;
			if (~isnan(cov_data.TRACER(i).cov(k,1)) && results.MET.qaqc.completeness(linenum) > 0.95 && results.MET.ust(linenum) > 0.3 && results.TRACER(i).qaqc.completeness(linenum) > 0.3 && results.TRACER(i).qaqc.SST(linenum) < 0.3)
                segment_idx = find(data_segments(1:end-1) <= results.time(linenum) & data_segments(2:end) > results.time(linenum), 1);
                cumul_cov(i,:,segment_idx) = cumul_cov(i,:,segment_idx) + cov_data.TRACER(i).cov(k,:);
                cumul_abs_cov(i,:,segment_idx) = cumul_abs_cov(i,:,segment_idx) + abs(cov_data.TRACER(i).cov(k,:));
			end
		end
		
	end

	progress_bar(length(cov_files), n);
	
end

% find lag times from cumulated covariances
fprintf('step 2 of 3: determining lagtimes\n');
for (i = 1:length(results.TRACER))

    if (~isempty(tracer_indices) && sum(i == tracer_indices) == 0)
        % skip this tracer
        continue;
    end

    for (s = 1:length(data_segments)-1)

        idx = find(results.time >= data_segments(s) & results.time < data_segments(s+1));
        lag_samples = find_covariance_peak(cumul_abs_cov(i,:,s), MAX_LAG, LAG_SEARCH_RANGE, COVPEAK_FILTER_LENGTH); % lag in samples
        lagtime = lag_samples/SAMPLING_RATE_SONIC; % seconds
        results.TRACER(i).lagtime2(idx) = lagtime;

        % save plots of cumulated covariances
        subplot(2,2,1); plot(-MAX_LAG:MAX_LAG, cumul_abs_cov(i,:,s)); ylabel('cumul. abs. cov.');
        title(sprintf('%s, segment %d: %.1f s', results.TRACER(i).name, s, lagtime));
        subplot(2,2,3); plot((-LAG_SEARCH_RANGE:LAG_SEARCH_RANGE), cumul_abs_cov(i,(MAX_LAG+1-LAG_SEARCH_RANGE):(MAX_LAG+1+LAG_SEARCH_RANGE),s)); xlabel('lag [samples]'); ylabel('cumul. abs. cov.');
        subplot(2,2,2); plot(-MAX_LAG:MAX_LAG, cumul_cov(i,:,s)); ylabel('cumul. covariance');
        subplot(2,2,4); plot((-LAG_SEARCH_RANGE:LAG_SEARCH_RANGE), cumul_cov(i,(MAX_LAG+1-LAG_SEARCH_RANGE):(MAX_LAG+1+LAG_SEARCH_RANGE),s)); xlabel('lag [samples]'); ylabel('cumul. covariance');
        output_file = sprintf('%s/figures/cumulated_cov%d_%d', output_folder, i, s);
        print(output_file, '-dpng');

    end
	progress_bar(length(results.TRACER), i);
end

% loop over all covariance files
fprintf('step 3 of 3: calculating fluxes\n');
for (n = 1:length(cov_files))

	% load covariance file
	cov_data = load([cov_files(n).path '/' cov_files(n).name]);
	
	% loop over all tracers
	for (i = 1:length(cov_data.TRACER))
		
        if (~isempty(tracer_indices) && sum(i == tracer_indices) == 0)
            % skip this tracer
            continue;
        end
		
		% find position in results file corresponding to current covariance file
		[~, results_idx] = min(abs(cov_files(n).date - results.time));
		
        % calculate flux using lagtimes of cumulated covariances
 		for (k = 1:intervals_per_day)
            linenum = results_idx-1+k;
            flux_conversion_factor = 1/(results.MET.T(linenum)/273.15*22.4e-3);
            if (~isnan(results.TRACER(i).lagtime2(linenum)))
                lag_samples = int32(results.TRACER(i).lagtime2(linenum)*SAMPLING_RATE_SONIC);
                results.TRACER(i).flux2(linenum) = cov_data.TRACER(i).cov(k, MAX_LAG+1-lag_samples)*flux_conversion_factor;
                % WPL correction
                if (APPLY_WPL_CORRECTION)
                    flux_individual = flux_individual + results.TRACER(i).wpl_corr(linenum);
                end
            end
        end
        
    end

	progress_bar(length(cov_files), n);

end

% save results to file
save([output_folder '/results.mat'], '-struct', 'results');
