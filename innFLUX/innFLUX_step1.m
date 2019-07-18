% innFLUX_step1.m

% the sonic data files are expected to cover one day of data each

results.pocedure_version = '2018.07.23';

% load parameters
innFLUX_parameters;

% store parameters
results.parameters.SONIC_ORIENTATION = SONIC_ORIENTATION;           % degrees
results.parameters.SENSOR_HEIGHT = SENSOR_HEIGHT;                   % sensor height in meters above roughness height
results.parameters.WINDOW_LENGTH = WINDOW_LENGTH;                   % samples
results.parameters.SAMPLING_RATE_SONIC = SAMPLING_RATE_SONIC;       % per second; sonic anemometer sampling rate
results.parameters.SAMPLING_RATE_TRACER = SAMPLING_RATE_TRACER;     % per second; gas analyzer sampling rate
results.parameters.DISJUNCT_EC = DISJUNCT_EC;                       % if 1, disjunct eddy covariance is applied
results.parameters.DETREND_TRACER_SIGNAL = DETREND_TRACER_SIGNAL;   % if 1, tracer signal detrending is applied
results.parameters.MAX_LAG = MAX_LAG;                               % samples; maximum lag for covariance function
results.parameters.LAG_SEARCH_RANGE = LAG_SEARCH_RANGE;             % samples; maximum lag for lag-time search
results.parameters.COVPEAK_FILTER_LENGTH = COVPEAK_FILTER_LENGTH;   % samples; length of filter applied to covariance peak prior to lag-time determination
results.parameters.NUM_FREQ_BINS = NUM_FREQ_BINS;                   % number of logarithmically spaced frequency bins for cospectra
results.parameters.APPLY_TILT_CORRECTION = APPLY_TILT_CORRECTION;   % if 1, sonic tilt correction is applied
results.parameters.APPLY_WPL_CORRECTION = APPLY_WPL_CORRECTION;     % if 1, WPL correction is applied
results.parameters.COMPLETENESS_THRESHOLD = COMPLETENESS_THRESHOLD; % threshold value for completeness of input data (0.0 - 1.0)

% store current time
proc_start_time = datenum(clock());
fprintf('procedure started at %s\n', datestr(proc_start_time));

% get list of input files
sonic_files_list = [];
for (n = 1:length(sonic_files_folders))
	sonic_files_list = [sonic_files_list; list_files_datestring(sonic_files_folders{n}, {'w'}, '.mat', 'yyyymmdd')];
end
tracer_files_list = [];
if (~isempty(tracer_files_folders))
    for (n = 1:length(tracer_files_folders))
        tracer_files_list = [tracer_files_list; list_files_datestring(tracer_files_folders{n}, {tracer_files_prefix}, '.mat', 'yyyymmdd')];
    end
else
    for (n = 1:length(tracer_files))
        [pathstr, name, ext] = fileparts(tracer_files{n});
        tracer_file.name = name;
        tracer_file.path = pathstr;
        tracer_files_list = [tracer_files_list; tracer_file];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sonic_files_list = sonic_files_list(134:178); % 0513-0626
%sonicdata =sonic_files_list = sonic_files_list(191:252); % 0710-0909
%sonic_files_list = sonic_files_list(182:293); % 0701-1020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('innFLUX_step1: found %d tracer files and %d sonic files\n', length(tracer_files_list), length(sonic_files_list));

% load tilt correction file
if (APPLY_TILT_CORRECTION)
    fprintf('using tilt correction file %s\n', tilt_correction_filepath);
    tilt_correction = load(tilt_correction_filepath);
end

% load pressure file
pressuredata = [];
if (~isempty(pressure_filepath))
    fprintf('using pressure file %s\n', pressure_filepath);
    pressuredata = load(pressure_filepath);
end

linenum = 0;
intervals_per_day = int32(floor(86400*SAMPLING_RATE_SONIC/WINDOW_LENGTH));
len = intervals_per_day*length(sonic_files_list);
upsampling_factor = SAMPLING_RATE_SONIC/SAMPLING_RATE_TRACER;

% load first tracer file
file_index = 1; % index of current tracer data file
if (~isempty(tracer_files_list))
    fprintf('tracer file %s\n', tracer_files_list(file_index).name);
    tracerdata_fullfile = load([tracer_files_list(file_index).path '/' tracer_files_list(file_index).name]);
    % upsample data to sonic sampling rate if necessary
    if (SAMPLING_RATE_TRACER < SAMPLING_RATE_SONIC && ~DISJUNCT_EC)
        tracerdata_fullfile.data = upfirdn(tracerdata_fullfile.data, ones(upsampling_factor, size(tracerdata_fullfile.data, 2)), upsampling_factor); % step-function
        tracerdata_fullfile.data(:,1) = (tracerdata_fullfile.data(1,1) + (0:(size(tracerdata_fullfile.data, 1)-1))/SAMPLING_RATE_SONIC/86400)'; % first column contains timestamp
    end
end

% check output folder
if (exist(output_folder, 'dir') == 0)
    fprintf('output folder not found: %s\n', output_folder);
    return;
end
% create folder for covariance files
if (exist([output_folder '/cov'], 'dir') == 7)
    % cov folder exists, delete
    rmdir([output_folder '/cov'], 's');
end
if (mkdir([output_folder '/cov']) == 0)
    fprintf('could not create folder for covariance files: %s\n', [output_folder '/cov']);
    return;
end

% prepare results structure
num_tracers = 0; % number of tracers for which to calculate fluxes
if (~isempty(tracer_files_list))
    num_tracers = size(tracerdata_fullfile.data, 2) - 1;
end
results.time = NaN(len,1);
results.hour = NaN(len,1); % hour of day
results.freq = []; % frequency axis of cospectra [1/s]
results.freq_scaled = []; % frequency axis of scaled cospectra
results.MET.uvw = NaN(len,3); % mean wind speed vector with u component rotated into direction of mean wind [m/s]
results.MET.std_uvw = NaN(len,3); % standard deviations for wind speed components u, v, w [m/s]
results.MET.hws = NaN(len,1); % mean horizontal wind speed [m/s]
results.MET.wdir = NaN(len,1); % mean horizontal wind direction [deg]
tilt.P = [];
results.MET.tilt = repmat(tilt,len,1); % wind tilt correction matrix
results.MET.uw = NaN(len,1); % covariance of along-wind and vertical wind component, <u'w'> [m^2/s^2]
results.MET.vw = NaN(len,1); % covariance of cross-wind and vertical wind component, <v'w'> [m^2/s^2]
results.MET.uv = NaN(len,1); % covariance of along-wind and cross-wind component, <u'v'> [m^2/s^2]
results.MET.uu = NaN(len,1); % auto-covariance of along-wind component, <u'u'> [m^2/s^2]
results.MET.vv = NaN(len,1); % auto-covariance of cross-wind component, <v'v'> [m^2/s^2]
results.MET.ww = NaN(len,1); % auto-covariance of vertical wind component, <w'w'> [m^2/s^2]
results.MET.ust = NaN(len,1); % friction velocity, u* [m/s]
results.MET.T = NaN(len,1); % mean temperature [K]
results.MET.std_T = NaN(len,1); % stddev of temperature [K]
results.MET.wT = NaN(len,1); % temperature flux, <w'T'> [K*m/s]
results.MET.L = NaN(len,1); % Obukhov length [m]
results.MET.zoL = NaN(len,1); % stability parameter, z/L
results.MET.cospec_wT = NaN(len,NUM_FREQ_BINS); % cospectra for wT (f*Co)
results.MET.cospec_wT_scaled = NaN(len,NUM_FREQ_BINS); % scaled cospectra for wT (f*Co/cov)
results.MET.p = NaN(len,1); % pressure [hPa]
results.MET.theta = NaN(len,1); % potential temperature [K]
results.MET.theta_v = NaN(len,1); % virtual potential temperature [K]
results.MET.wtheta = NaN(len,1); % potential temperature (heat) flux, <w'theta'> [K*m/s]
results.MET.wtheta_v = NaN(len,1); % vitual potential temperature (buoyancy) flux, <w'theta_v'> [K*m/s]
results.MET.qaqc.completeness = NaN(len,1); % fraction of sonic data used in this interval
results.MET.qaqc.SST_wT = NaN(len,1); % steady state test for wT
results.MET.qaqc.ITC_w = NaN(len,1); % integral turbulence characteristics test for w
results.MET.qaqc.ITC_u = NaN(len,1); % integral turbulence characteristics test for u
results.MET.qaqc.ITC_T = NaN(len,1); % integral turbulence characteristics test for T
results.MET.qaqc.cospec_wT_integral = NaN(len,1); % integrated cospectrum for wT
if (~isempty(irga_columns))
    IRGA.name = [];
    IRGA.mean = NaN(len,1); % mean concentration [ppbv]
    IRGA.std = NaN(len,1); % stddev of concentration
    IRGA.flux = NaN(len,1); % flux [nmol/m^2/s]
    IRGA.cospec = NaN(len,NUM_FREQ_BINS); % cospectrum: f*Co
    IRGA.cospec_scaled = NaN(len,NUM_FREQ_BINS); % scaled cospectrum f*Co/cov
    IRGA.qaqc.completeness = NaN(len,1); % fraction of tracer data used in this interval
    IRGA.qaqc.SST = NaN(len,1); % steady state test for tracer
    IRGA.qaqc.flux_SNR = NaN(len,1); % flux signal to noise ratio
    IRGA.qaqc.flux_noise_std = NaN(len,1); % stddev of flux noise far off integral time scale
    IRGA.qaqc.flux_noise_mean = NaN(len,1); % mean flux noise far off integral time scale
    IRGA.qaqc.flux_noise_rmse = NaN(len,1); % RMSE of flux noise far off integral time scale
    IRGA.qaqc.random_error_FS = NaN(len,1); % random error as described by Finkelstein and Sims (2001)
    IRGA.qaqc.random_error_noise = NaN(len,1); % random error noise estimated according to Mauder 2013
    IRGA.qaqc.random_flux = NaN(len,1); % random flux level estimated by random shuffle criterium (Billesbach 2011)
    IRGA.qaqc.cospec_integral = NaN(len,1); % integrated cospectrum
    results.IRGA = repmat(IRGA, 1, length(irga_columns)); % struct containing results for each tracer
    clear IRGA;
    for (i = 1:length(irga_columns))
        results.IRGA(i).name = irga_names{i};
    end
end
if (~isempty(tracer_files_list))
    TRACER.name = []; % tracer name
    TRACER.mean = NaN(len,1); % mean tracer concentration [ppb]
    TRACER.std = NaN(len,1); % stddev of detrended tracer concentration [ppb]
    TRACER.lagtime1 = NaN(len,1); % [s]
    TRACER.flux1 = NaN(len,1); % flux from individual lagtime [nmol/m^2/s]
    TRACER.wpl_corr = NaN(len,1); % WPL (Webb) correction term [nmol/m^2/s]
    TRACER.cospec = NaN(len,NUM_FREQ_BINS); % cospectrum: f*Co
    TRACER.cospec_scaled = NaN(len,NUM_FREQ_BINS); % scaled cospectrum: f*Co/cov
    TRACER.qaqc.completeness = NaN(len,1); % fraction of tracer data used in this interval
    TRACER.qaqc.SST = NaN(len,1); % steady state test for tracer
    TRACER.qaqc.flux_SNR = NaN(len,1); % flux signal to noise ratio
    TRACER.qaqc.flux_noise_std = NaN(len,1); % stddev of flux noise far off integral time scale
    TRACER.qaqc.flux_noise_mean = NaN(len,1); % mean flux noise far off integral time scale
    TRACER.qaqc.flux_noise_rmse = NaN(len,1); % RMSE of flux noise far off integral time scale
    TRACER.qaqc.random_error_FS = NaN(len,1); % random error as described by Finkelstein and Sims (2001)
    TRACER.qaqc.random_error_noise = NaN(len,1); % random error noise estimated according to Mauder 2013
    TRACER.qaqc.random_flux = NaN(len,1); % random flux level estimated by random shuffle criterium (Billesbach 2011)
    TRACER.qaqc.cospec_integral = NaN(len,1); % integrated cospectrum
    results.TRACER = repmat(TRACER, 1, num_tracers); % struct containing results for each tracer
    clear TRACER;
    for (i = 1:num_tracers)
        results.TRACER(i).name = tracerdata_fullfile.header{i+1};
    end
end

% frequency axis for cospectra
f = (0:ceil(WINDOW_LENGTH/2))/WINDOW_LENGTH*SAMPLING_RATE_SONIC; % linear frequency
[results.freq, ~] = logBinSpectrum(f, zeros(size(f)), NUM_FREQ_BINS, f(2), f(end)); % log-spaced freuqency bins
[results.freq_scaled, ~] = logBinSpectrum(f, zeros(size(f)), NUM_FREQ_BINS, FREQ_BIN_MIN, FREQ_BIN_MAX); % log-spaced freuqency bins for scaled cospectra

for (n = 1:length(sonic_files_list))

	% covariance data to be stored to file
	for (i = 1:num_tracers)
		cov_data.TRACER(i).cov = NaN(intervals_per_day, 2*MAX_LAG+1);
	end

	% load sonic file
	fprintf('sonic file %s\n', sonic_files_list(n).name);
	sonicdata_fullday = load([sonic_files_list(n).path '/' sonic_files_list(n).name]); % columns: time, ...

	% loop over all intervals (of WINDOW_LENGTH) of the day
	for (day_interval = 1:intervals_per_day)

		% define interval
		interval_start_time = floor(sonicdata_fullday.data(1,1)) + double(day_interval-1)*WINDOW_LENGTH/SAMPLING_RATE_SONIC/86400;
		interval_end_time = interval_start_time + (WINDOW_LENGTH-1)/SAMPLING_RATE_SONIC/86400;

        % select sonic data interval
		[~, start_idx] = min(abs(sonicdata_fullday.data(:,1) - interval_start_time));
		[~, end_idx] = min(abs(sonicdata_fullday.data(:,1) - interval_end_time));
		sonicdata = sonicdata_fullday.data(start_idx:end_idx, :);
        
		% store timestamp
		linenum = linenum + 1;
		results.time(linenum) = interval_start_time;
		results.hour(linenum) = 24*(interval_start_time - floor(interval_start_time));
        
        % check sonic flag and set flagged data to NaN
        flag_sonic = sonicdata(:,7) ~= 0;
        sonicdata(flag_sonic, 3:6) = NaN;

        % check completeness of sonic data
        completeness_sonicdata = sum(isfinite(sonicdata(:,3)))/WINDOW_LENGTH;

        if (completeness_sonicdata >= COMPLETENESS_THRESHOLD)
        
            % sonic tilt correction
            mean_uvw = nanmean(sonicdata(:, 3:5), 1);
            mean_wind_dir = mod( 180/pi*atan2( -mean_uvw(2), mean_uvw(1) ) + SONIC_ORIENTATION, 360 );
            if (APPLY_TILT_CORRECTION)
                P_tilt = tilt_correction.data(mod(int32(mean_wind_dir),360) + 1).P;
                wind_detilted = (P_tilt*[sonicdata(:,3)'; sonicdata(:,4)'; sonicdata(:,5)'])'; % untilted u, v and w
            else
                P_tilt = [];
                wind_detilted = [sonicdata(:,3) sonicdata(:,4) sonicdata(:,5)]; % untilted u, v and w
            end

            % mean horizontal wind speed and direction
            % note: wind direction is the angle where the wind is coming from, with north = 0 and east = 90 degrees
            % sonic anemometer coordinate system: see CSAT3 manual
            mean_uvw = nanmean(wind_detilted, 1);
            mean_wind_speed = sqrt(mean_uvw(1)^2 + mean_uvw(2)^2 + mean_uvw(3)^2);
            mean_wind_dir = mod( 180/pi*atan2( -mean_uvw(2), mean_uvw(1) ) + SONIC_ORIENTATION, 360 );

            % rotate wind vectors so that u component points along mean wind
            angle = atan2(mean_uvw(2), mean_uvw(1));
            R = [cos(-angle), -sin(-angle); sin(-angle), cos(-angle)];
            wind_detilted(:,1:2) = (R*(wind_detilted(:,1:2)'))';
            % recalculate mean wind
            mean_uvw = nanmean(wind_detilted, 1);
            
            % interpolate wind and scalar data (keep NaN and extrapolate using NaN)
            interp_time = interval_start_time + (0:(WINDOW_LENGTH-1))'/SAMPLING_RATE_SONIC/3600/24;
            interp_u = interp1(sonicdata(:,1), wind_detilted(:,1), interp_time, 'nearest', NaN);
            interp_v = interp1(sonicdata(:,1), wind_detilted(:,2), interp_time, 'nearest', NaN);
            interp_w = interp1(sonicdata(:,1), wind_detilted(:,3), interp_time, 'nearest', NaN);
            interp_Tsonic = interp1(sonicdata(:,1), sonicdata(:,6), interp_time, 'nearest', NaN) + 273.15;
            
            % check data completeness
            completeness_T = sum(isfinite(interp_Tsonic))/WINDOW_LENGTH;

            % Reynolds decomposition of wind with detrending (e.g. u = u_mean + u_prime)
            u_prime = nandetrend(interp_u);
            v_prime = nandetrend(interp_v);
            w_prime = nandetrend(interp_w);

            % standard deviation of wind
            std_uvw = [nanstd(u_prime) nanstd(v_prime) nanstd(w_prime)];

            % replace NaN by mean value (0)
            u_prime = replace_nan(u_prime, 0);
            v_prime = replace_nan(v_prime, 0);
            w_prime = replace_nan(w_prime, 0);
            
            % calculate covariances of wind components
            uw = xcov(u_prime, w_prime, 0, 'unbiased');
            vw = xcov(v_prime, w_prime, 0, 'unbiased');
            uv = xcov(u_prime, v_prime, 0, 'unbiased');
            uu = xcov(u_prime, u_prime, 0, 'unbiased');
            vv = xcov(v_prime, v_prime, 0, 'unbiased');
            ww = xcov(w_prime, w_prime, 0, 'unbiased');
            
            % calculate friction velocity, u*
            u_star = (uw^2 + vw^2)^0.25;
                        
            % calculate temperature T from virtual sonic temperature T_sonic and water vapor concentration
            % T_sonic = T*(1 + 0.32*w), were w is the volume mixing ratio of water vapor in air
            if (~isempty(irga_columns) && irga_H2O_index > 0)
                
                % check IRGA flag and set flagged data to NaN
                flag_irga = sonicdata(:, irga_flag_column) ~= 0;
                sonicdata(flag_irga, irga_columns) = NaN;

                interp_H2O = interp1(sonicdata(:,1), sonicdata(:,irga_columns(irga_H2O_index)), interp_time, 'nearest', NaN);
                completeness_H2O = sum(isfinite(interp_H2O))/WINDOW_LENGTH;
                if (completeness_T >= COMPLETENESS_THRESHOLD && completeness_H2O >= COMPLETENESS_THRESHOLD)
                    interp_T = interp_Tsonic./(1 + 0.32*1e-9*interp_H2O);
                    T_mean = nanmean(interp_T);
                    T_prime = nandetrend(interp_T);
                    std_T = nanstd(T_prime);
                    T_prime = replace_nan(T_prime, 0);
                end
                
            else
                
                % H2O concentration not available, use T_sonic
                completeness_H2O = 0;
                interp_T = interp_Tsonic;
                T_mean = nanmean(interp_Tsonic);
                T_prime = nandetrend(interp_Tsonic);
                std_T = nanstd(interp_Tsonic);
                T_prime = replace_nan(T_prime, 0);
                
            end
            
            % IRGA tracers
            if (~isempty(irga_columns))
                
                % check IRGA flag and set flagged data to NaN
                flag_irga = sonicdata(:, irga_flag_column) ~= 0;
                sonicdata(flag_irga, irga_columns) = NaN;
                
                for (i = 1:length(irga_columns))
                
                    interp_c = interp1(sonicdata(:,1), sonicdata(:,irga_columns(i)), interp_time, 'nearest', NaN);
                    completeness_IRGA = sum(isfinite(interp_c))/WINDOW_LENGTH;
                    
                    if (completeness_IRGA >= COMPLETENESS_THRESHOLD)
                        
                        c_mean = nanmean(interp_c);
                        c_prime = interp_c - c_mean;
                        std_c = nanstd(c_prime);
                        c_prime = replace_nan(c_prime, 0);
                        
                        % calculate cross-covariance
                        cov_wc = xcov(w_prime, c_prime, MAX_LAG, 'unbiased');

                        % calculate flux in nmol/m^2/s
                        flux_conversion_factor = 1/(T_mean/273.15*22.4e-3);
                        flux = cov_wc(MAX_LAG + 1)*flux_conversion_factor;

                        % calculate cospectra
                        cospec = real( conj(fft(w_prime)) .* fft(c_prime) )';
                        cospec = cospec/(WINDOW_LENGTH^2);
                        cospec_wc_integral = sum(cospec);
                        [~, cospec_wc] = logBinSpectrum(f, f.*cospec(1:length(f)), NUM_FREQ_BINS, f(2), f(end)); % log-spaced freuqency bins
                        [~, cospec_wc_scaled] = logBinSpectrum(f*SENSOR_HEIGHT/mean_wind_speed, f.*cospec(1:length(f))/cospec_wc_integral, NUM_FREQ_BINS, FREQ_BIN_MIN, FREQ_BIN_MAX);

                        % steady state test according to Foken and Wichura (1996)
                        steady_state_test_wc = steady_state_test(w_prime, c_prime);

                        % flux detection limit: flux noise criterium using STD noise of covariance between +/- 160-180s (Spirig et al. 2005)
                        idx = [(MAX_LAG + 1 - 180*SAMPLING_RATE_SONIC):(MAX_LAG + 1 - 160*SAMPLING_RATE_SONIC), (MAX_LAG + 1 + 160*SAMPLING_RATE_SONIC):(MAX_LAG + 1 + 180*SAMPLING_RATE_SONIC)];
                        flux_noise_mean = nanmean(cov_wc(idx))*flux_conversion_factor;
                        flux_noise_std = nanstd(cov_wc(idx))*flux_conversion_factor;

                        % flux detection limit: flux noise criterium using RMSE noise of covariance between +/- 160-180s
                        idx_left = (MAX_LAG + 1 - 180*SAMPLING_RATE_SONIC):(MAX_LAG + 1 - 160*SAMPLING_RATE_SONIC);
                        idx_right = (MAX_LAG + 1 + 160*SAMPLING_RATE_SONIC):(MAX_LAG + 1 + 180*SAMPLING_RATE_SONIC);
                        flux_noise_rmse = sqrt( 0.5*( nanstd(cov_wc(idx_left))^2 + nanmean(cov_wc(idx_left))^2 + nanstd(cov_wc(idx_right))^2 + nanmean(cov_wc(idx_right))^2 ) )*flux_conversion_factor;
                        
                        % estimate random error as described by Finkelstein and Sims (2001)
                        random_error_FS = sqrt( 1/WINDOW_LENGTH * ( sum( xcov(w_prime, w_prime, 100, 'unbiased').*xcov(c_prime, c_prime, 100, 'unbiased') ) + sum( xcov(w_prime, c_prime, 100, 'unbiased').*xcov(c_prime, w_prime, 100, 'unbiased') ) ) )*flux_conversion_factor;
                        
                        % flux detection limit: random shuffle (Billesbach 2011)
                        for (irand = 1:10)
                            w_rand = w_prime(randperm(length(w_prime)));
                            xcov_rand(irand) = xcov(w_rand, c_prime, 0, 'unbiased');
                        end
                        random_flux = nanstd(xcov_rand)*flux_conversion_factor;

                        % estimate white noise using autocovariance (Mauder et al. 2013???)
                        % TODO: linear fit instead of interp1
                        autocov_c = xcov(c_prime, c_prime, 5, 'unbiased');
                        white_noise_c = ( sqrt(autocov_c(6)) - sqrt(interp1(7:11, autocov_c(7:11), 6, 'linear', 'extrap')) )*flux_conversion_factor;
                        random_error_noise = white_noise_c * nanstd(w_prime) / sqrt(WINDOW_LENGTH);

                        % store results for individual IRGA tracer
                        results.IRGA(i).mean(linenum) = c_mean;
                        results.IRGA(i).std(linenum) = std_c;
                        results.IRGA(i).flux(linenum) = flux;
                        results.IRGA(i).cospec(linenum,:) = cospec_wc;
                        results.IRGA(i).cospec_scaled(linenum,:) = cospec_wc_scaled;
                        results.IRGA(i).qaqc.completeness(linenum) = completeness_IRGA;
                        results.IRGA(i).qaqc.SST(linenum) = steady_state_test_wc;
                        results.IRGA(i).qaqc.flux_SNR(linenum) = abs(flux - flux_noise_mean)/flux_noise_std;
                        results.IRGA(i).qaqc.flux_noise_mean(linenum) = flux_noise_mean;
                        results.IRGA(i).qaqc.flux_noise_std(linenum) = flux_noise_std;
                        results.IRGA(i).qaqc.flux_noise_rmse(linenum) = flux_noise_rmse;
                        results.IRGA(i).qaqc.random_error_FS(linenum) = random_error_FS;
                        results.IRGA(i).qaqc.random_flux(linenum) = random_flux;
                        results.IRGA(i).qaqc.random_error_noise(linenum) = random_error_noise;
                        results.IRGA(i).qaqc.cospec_integral(linenum) = cospec_wc_integral;
                    end
                    
                end
                            
            end

            % store H2O concentration and flux for later
            if (~isempty(irga_columns) && irga_H2O_index > 0)
                H2O_mean = results.IRGA(irga_H2O_index).mean(linenum);
                flux_H2O = results.IRGA(irga_H2O_index).flux(linenum);
            else
                H2O_mean = NaN;
                flux_H2O = NaN;
            end

            if (completeness_T >= COMPLETENESS_THRESHOLD)
                
                % calculate buoyancy flux <w'T'>
                wT = xcov(w_prime, T_prime, 0, 'unbiased');

                % calculate cospectra for wT
                cospec = real( conj(fft(w_prime)) .* fft(T_prime) )';
                cospec = cospec/(WINDOW_LENGTH^2);
                cospec_wT_integral = sum(cospec);
                [~, cospec_wT] = logBinSpectrum(f, f.*cospec(1:length(f)), NUM_FREQ_BINS, f(2), f(end)); % log-spaced freuqency bins
                [~, cospec_wT_scaled] = logBinSpectrum(f*SENSOR_HEIGHT/mean_wind_speed, f.*cospec(1:length(f))/cospec_wT_integral, NUM_FREQ_BINS, FREQ_BIN_MIN, FREQ_BIN_MAX);

                % steady state test for wT according to Foken and Wichura (1996)
                steady_state_test_wT = steady_state_test(w_prime, T_prime);
                
            else
                wT = NaN;
                cospec_wT = NaN(1, NUM_FREQ_BINS);
                cospec_wT_scaled = NaN(1, NUM_FREQ_BINS);
                cospec_wT_integral = NaN;
                steady_state_test_wT = NaN;
            end

            % calculate Obukhov length, L
            L = -T_mean*u_star^3/(0.4*9.81*wT);
            
            % calculate stability parameter, z/L
            zoL = SENSOR_HEIGHT/L; % z over L

            % calculate pressure-dependent values
            if (isempty(pressuredata))
                p_mean = DEFAULT_PRESSURE;
            else
                start_idx = find(pressuredata.time >= interval_start_time, 1);
                end_idx = find(pressuredata.time >= interval_end_time, 1);
                p_mean = nanmean(pressuredata.p(start_idx:end_idx));
            end
            p_0 = 1000; % reference pressure for potential temperature [hPa]
            if ((completeness_T >= COMPLETENESS_THRESHOLD) && isfinite(p_mean))
                theta_mean = T_mean*(p_0/p_mean)^(0.286);
                theta_prime = T_prime*(p_0/p_mean)^(0.286);
                wtheta = xcov(w_prime, theta_prime, 0, 'unbiased');
                if (completeness_H2O >= COMPLETENESS_THRESHOLD)
                    theta_v = interp_T*(p_0/p_mean)^(0.286).*(1 + 0.38*1e-9*interp_H2O);
                    theta_v_mean = nanmean(theta_v);
                    theta_v_prime = nandetrend(theta_v);
                    theta_v_prime = replace_nan(theta_v_prime, 0);
                    wtheta_v = xcov(w_prime, theta_v_prime, 0, 'unbiased');
                else
                    theta_v_mean = NaN;
                    wtheta_v = NaN;
                end
            else
                theta_mean = NaN;
                theta_v_mean = NaN;
                wtheta = NaN;
                wtheta_v = NaN;
            end
            
            % test on developed turbulent conditions
            % Foken and Wichura (1996)
            % calculate integral turbulence characteristics based on flux-variance similarity
            if (zoL < -1)
                norm_sigma_w_model = 2*(-zoL)^(1/6);
                norm_sigma_u_model = 2.83*(-zoL)^(1/6);
                norm_sigma_T_model = (-zoL)^(-1/3);
            elseif (zoL < -0.0625)
                norm_sigma_w_model = 2*(-zoL)^(1/8);
                norm_sigma_u_model = 2.83*(-zoL)^(1/8);
                norm_sigma_T_model = (-zoL)^(-1/4);
            elseif (zoL < 0)
                norm_sigma_w_model = 1.41;
                norm_sigma_u_model = 1.99;
                norm_sigma_T_model = 0.5*(-zoL)^(-1/2);
            else
                norm_sigma_w_model = NaN;
                norm_sigma_u_model = NaN;
                norm_sigma_T_model = NaN;
            end
            T_star = -wT/u_star;
            norm_sigma_w_meas = nanstd(w_prime)/u_star;
            norm_sigma_u_meas = nanstd(u_prime)/u_star;
            norm_sigma_T_meas = nanstd(T_prime)/T_star;
            ITC_w = abs((norm_sigma_w_model - norm_sigma_w_meas)/norm_sigma_w_model);
            ITC_u = abs((norm_sigma_u_model - norm_sigma_u_meas)/norm_sigma_u_model);
            ITC_T = abs((norm_sigma_T_model - norm_sigma_T_meas)/norm_sigma_T_model);

            % store results for sonic data
            results.MET.uvw(linenum,:) = mean_uvw;
            results.MET.std_uvw(linenum,:) = std_uvw;
            results.MET.hws(linenum) = mean_wind_speed;
            results.MET.wdir(linenum) = mean_wind_dir;
            tilt.P = P_tilt;
            results.MET.tilt(linenum) = tilt;
            results.MET.uw(linenum) = uw;
            results.MET.vw(linenum) = vw;
            results.MET.uv(linenum) = uv;
            results.MET.uu(linenum) = uu;
            results.MET.vv(linenum) = vv;
            results.MET.ww(linenum) = ww;
            results.MET.ust(linenum) = u_star;
            results.MET.T(linenum) = T_mean;
            results.MET.std_T(linenum) = std_T;
            results.MET.wT(linenum) = wT;
            results.MET.L(linenum) = L;
            results.MET.zoL(linenum) = zoL;
            results.MET.cospec_wT(linenum,:) = cospec_wT;
            results.MET.cospec_wT_scaled(linenum,:) = cospec_wT_scaled;
            if (~isempty(pressuredata))
                results.MET.p(linenum) = p_mean;
                results.MET.theta(linenum) = theta_mean;
                results.MET.theta_v(linenum) = theta_v_mean;
                results.MET.wtheta(linenum) = wtheta;
                results.MET.wtheta_v(linenum) = wtheta_v;
            end
            results.MET.qaqc.completeness(linenum) = completeness_sonicdata;
            results.MET.qaqc.SST_wT(linenum) = steady_state_test_wT;
            results.MET.qaqc.ITC_w(linenum) = ITC_w;
            results.MET.qaqc.ITC_u(linenum) = ITC_u;
            results.MET.qaqc.ITC_T(linenum) = ITC_T;
            results.MET.qaqc.cospec_wT_integral(linenum) = cospec_wT_integral;

            % process tracer data
            if (~isempty(tracer_files_list))
            
                % load next tracer file if necessary
                while (file_index < length(tracer_files_list) && tracerdata_fullfile.data(end,1) < interval_end_time - 0.1*WINDOW_LENGTH/SAMPLING_RATE_SONIC/86400)
                    % load next tracer file
                    file_index = file_index + 1;
                    fprintf('tracer file %s\n', tracer_files_list(file_index).name);
                    tracerdata_fullfile = load([tracer_files_list(file_index).path '/' tracer_files_list(file_index).name]);
                    % upsample data to sonic sampling rate if necessary
                    if (SAMPLING_RATE_TRACER < SAMPLING_RATE_SONIC && ~DISJUNCT_EC)
                        tracerdata_fullfile.data = upfirdn(tracerdata_fullfile.data, ones(upsampling_factor, size(tracerdata_fullfile.data, 2)), upsampling_factor); % step-function
                        tracerdata_fullfile.data(:,1) = (tracerdata_fullfile.data(1,1) + (0:(size(tracerdata_fullfile.data, 1)-1))/SAMPLING_RATE_SONIC/86400)'; % first column contains timestamp
                    end
                end

                % select data interval
                start_idx = find(tracerdata_fullfile.data(:,1) >= interval_start_time, 1);
                end_idx = find(tracerdata_fullfile.data(:,1) >= interval_end_time, 1);
                tracerdata.time = tracerdata_fullfile.data(start_idx:end_idx, 1);
                tracerdata.ppb = tracerdata_fullfile.data(start_idx:end_idx, 2:end);

                % check completeness of tracer data
                if (DISJUNCT_EC)
                    completeness_tracerdata = sum(isfinite(tracerdata.ppb),1)/WINDOW_LENGTH*upsampling_factor;
                else
                    completeness_tracerdata = sum(isfinite(tracerdata.ppb),1)/WINDOW_LENGTH;
                end

                % iterate over all tracers
                for (i = 1:num_tracers)

                    if (~isempty(tracer_indices) && sum(i == tracer_indices) == 0)
                        % skip this tracer
                        continue;
                    end

                    if(completeness_tracerdata(i) >= COMPLETENESS_THRESHOLD)

                        % standard deviation
                        std_c = nanstd(tracerdata.ppb(:, i));

                        if (DISJUNCT_EC)

                            % detrend / subtract mean and upsample by gap-filling with zeros
                            if (DETREND_TRACER_SIGNAL)
                                upsampled_c.ppb = upsample(nandetrend(tracerdata.ppb(:, i)), upsampling_factor);
                            else
                                upsampled_c.ppb = upsample(tracerdata.ppb(:, i) - nanmean(tracerdata.ppb(:, i), 1), upsampling_factor);
                            end
                            upsampled_c.time = (tracerdata.time(1,1) + (0:(size(upsampled_c.ppb, 1)-1))/SAMPLING_RATE_SONIC/86400)';

                            % interpolate concentration data (keep NaN and extrapolate using NaN)
                            interp_c = interp1(upsampled_c.time, upsampled_c.ppb, interp_time, 'nearest', NaN);

                            % Reynolds decomposition (note: mean already subtracted from interp_c)
                            c_mean = nanmean(tracerdata.ppb(:, i), 1);
                            c_prime = interp_c;

                            % replace NaN by mean value (0)
                            c_prime = replace_nan(c_prime, 0);

                            % calculate cross-covariance
                            cov_wc = xcov(w_prime, c_prime, MAX_LAG, 'unbiased')*upsampling_factor;

                        else

                            % interpolate concentration data (keep NaN and extrapolate using NaN)
                            interp_c = interp1(tracerdata.time, tracerdata.ppb(:, i), interp_time, 'nearest', NaN);

                            % Reynolds decomposition
                            c_mean = nanmean(interp_c, 1);
                            if (DETREND_TRACER_SIGNAL)
                                c_prime = nandetrend(interp_c);
                            else
                                c_prime = interp_c - c_mean;
                            end

                            % replace NaN by mean value (0)
                            c_prime = replace_nan(c_prime, 0);

                            % calculate cross-covariance
                            cov_wc = xcov(w_prime, c_prime, MAX_LAG, 'unbiased');

                        end

                        % find lag time from covariance peak
                        lag_samples = find_covariance_peak(cov_wc, MAX_LAG, LAG_SEARCH_RANGE, COVPEAK_FILTER_LENGTH); % lag in samples
                        lagtime_individual = lag_samples/SAMPLING_RATE_SONIC; % seconds

                        % calculate flux in nmol/m^2/s
                        flux_conversion_factor = 1/(T_mean/273.15*22.4e-3);
                        flux_individual = cov_wc(MAX_LAG + 1 - lag_samples)*flux_conversion_factor;

                        % WPL correction
                        wpl_correction = 1e-9*c_mean/(1 - 1e-9*H2O_mean)*flux_H2O; % nmol/m^2/s
                        if (APPLY_WPL_CORRECTION)
                            flux_individual = flux_individual + wpl_correction;
                        end

                        % calculate cospectra
                        if (~DISJUNCT_EC)
                            c_prime_shifted = circshift(c_prime, -lag_samples);
                            if (lag_samples > 0)
                                c_prime_shifted((end-lag_samples+1):end) = 0; 
                            elseif (lag_samples < 0)
                                c_prime_shifted(1:(-lag_samples)) = 0;
                            end
                            cospec = real( conj(fft(w_prime)) .* fft(c_prime_shifted) )';
                            cospec = cospec/(WINDOW_LENGTH^2);
                            cospec_wc_integral = sum(cospec);
                            [~, cospec_wc] = logBinSpectrum(f, f.*cospec(1:length(f)), NUM_FREQ_BINS, f(2), f(end)); % log-spaced freuqency bins
                            [~, cospec_wc_scaled] = logBinSpectrum(f*SENSOR_HEIGHT/mean_wind_speed, f.*cospec(1:length(f))/cospec_wc_integral, NUM_FREQ_BINS, FREQ_BIN_MIN, FREQ_BIN_MAX);
                        else
                            cospec_wc_integral = NaN;
                            cospec_wc = NaN(1, NUM_FREQ_BINS);
                            cospec_wc_scaled = NaN(1, NUM_FREQ_BINS);
                        end

                        % steady state test according to Foken and Wichura (1996)
                        steady_state_test_wc = steady_state_test(w_prime, c_prime);

                        % flux detection limit: flux noise criterium using STD noise of covariance between +/- 160-180s (Spirig et al. 2005)
                        idx = [(MAX_LAG + 1 - 180*SAMPLING_RATE_SONIC):(MAX_LAG + 1 - 160*SAMPLING_RATE_SONIC), (MAX_LAG + 1 + 160*SAMPLING_RATE_SONIC):(MAX_LAG + 1 + 180*SAMPLING_RATE_SONIC)];
                        flux_noise_mean = nanmean(cov_wc(idx))*flux_conversion_factor;
                        flux_noise_std = nanstd(cov_wc(idx))*flux_conversion_factor;

                        % flux detection limit: flux noise criterium using RMSE noise of covariance between +/- 160-180s
                        idx_left = (MAX_LAG + 1 - 180*SAMPLING_RATE_SONIC):(MAX_LAG + 1 - 160*SAMPLING_RATE_SONIC);
                        idx_right = (MAX_LAG + 1 + 160*SAMPLING_RATE_SONIC):(MAX_LAG + 1 + 180*SAMPLING_RATE_SONIC);
                        flux_noise_rmse = sqrt( 0.5*( nanstd(cov_wc(idx_left))^2 + nanmean(cov_wc(idx_left))^2 + nanstd(cov_wc(idx_right))^2 + nanmean(cov_wc(idx_right))^2 ) )*flux_conversion_factor;

                        % estimate random error as described by Finkelstein and Sims (2001)
                        random_error_FS = sqrt( 1/WINDOW_LENGTH * ( sum( xcov(w_prime, w_prime, 100, 'unbiased').*xcov(c_prime, c_prime, 100, 'unbiased') ) + sum( xcov(w_prime, c_prime, 100, 'unbiased').*xcov(c_prime, w_prime, 100, 'unbiased') ) ) )*flux_conversion_factor;

                        % flux detection limit: random shuffle (Billesbach 2011)
                        for (irand = 1:10)
                            w_rand = w_prime(randperm(length(w_prime)));
                            xcov_rand(irand) = xcov(w_rand, c_prime, 0, 'unbiased');
                        end
                        random_flux = nanstd(xcov_rand)*flux_conversion_factor;

                        % estimate white noise using autocovariance (Mauder et al. 2013???)
                        % TODO: linear fit instead of interp1
                        autocov_c = xcov(c_prime, c_prime, 5, 'unbiased');
                        white_noise_c = ( sqrt(autocov_c(6)) - sqrt(interp1(7:11, autocov_c(7:11), 6, 'linear', 'extrap')) )*flux_conversion_factor;
                        random_error_noise = white_noise_c * nanstd(w_prime) / sqrt(WINDOW_LENGTH);

                        % store results for individual tracer
                        results.TRACER(i).mean(linenum) = c_mean;
                        results.TRACER(i).std(linenum) = std_c;
                        results.TRACER(i).lagtime1(linenum) = lagtime_individual;
                        results.TRACER(i).flux1(linenum) = flux_individual;
                        results.TRACER(i).wpl_corr(linenum) = wpl_correction;
                        results.TRACER(i).cospec(linenum,:) = cospec_wc;
                        results.TRACER(i).cospec_scaled(linenum,:) = cospec_wc_scaled;
                        results.TRACER(i).qaqc.completeness(linenum) = completeness_tracerdata(i);
                        results.TRACER(i).qaqc.SST(linenum) = steady_state_test_wc;
                        results.TRACER(i).qaqc.flux_SNR(linenum) = abs(flux_individual - flux_noise_mean)/flux_noise_std;
                        results.TRACER(i).qaqc.flux_noise_mean(linenum) = flux_noise_mean;
                        results.TRACER(i).qaqc.flux_noise_std(linenum) = flux_noise_std;
                        results.TRACER(i).qaqc.flux_noise_rmse(linenum) = flux_noise_rmse;
                        results.TRACER(i).qaqc.random_error_FS(linenum) = random_error_FS;
                        results.TRACER(i).qaqc.random_flux(linenum) = random_flux;
                        results.TRACER(i).qaqc.random_error_noise(linenum) = random_error_noise;
                        results.TRACER(i).qaqc.cospec_integral(linenum) = cospec_wc_integral;
                        cov_data.TRACER(i).cov(day_interval,:) = cov_wc;

                    end

                end
                
            end

        end

    end
	
	% save covariances for whole day to file
    if (~isempty(tracer_files_list))
    	dv = datevec(sonicdata_fullday.data(1,1));
        cov_filename = sprintf('cov%d.%02d.%02d.mat', dv(1), dv(2), dv(3));
        save([output_folder '/cov/' cov_filename], '-struct', 'cov_data');
    end
		
end

% save results to file
results.file_creation_time = datestr(now);
save([output_folder '/results.mat'], '-struct', 'results');

% display time needed
proc_end_time = datenum(clock());
fprintf('finished at %s. procedure took %s\n', datestr(proc_end_time), datestr(proc_end_time - proc_start_time, 'HH:MM:SS'));
