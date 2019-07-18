% calc_angle_dependent_tilt_correction.m
% calculates for a large wind data set the tilt angle of a sonic anemometer relative to the mean streamlines for sectors of mean wind direction

sonic_files_folders{1} = 'Z:\procdata\617m_CPEC\2018-03\dayfiles\fieldversions\vA';
sonic_files_folders{2} = 'Z:\procdata\617m_CPEC\2018-04\dayfiles\fieldversions\vA';
sonic_files_folders{3} = 'Z:\procdata\617m_CPEC\2018-05\dayfiles\fieldversions\vA';
sonic_files_folders{4} = 'Z:\procdata\617m_CPEC\2018-06\dayfiles\fieldversions\vA';
% sonic_files_folders{5} = 'Z:\procdata\617m_CPEC\IOP2018_Summer\dayfiles\vD\dayfiles';
% sonic_files_folders{6} = 'Z:\procdata\617m_CPEC\IOP2018_Summer\dayfiles\vE\dayfiles';
% sonic_files_folders{7} = 'Z:\procdata\617m_CPEC\IOP2018_Summer\dayfiles\vF\dayfiles';
% sonic_files_folders{8} = 'Z:\procdata\617m_CPEC\IOP2018_Summer\dayfiles\vG\dayfiles';
% sonic_files_folders{9} = 'Z:\procdata\617m_CPEC\IOP2018_Summer\dayfiles\vH\dayfiles';

% for G:\CPEC_preproc_2019-03-19 add following data
sonic_files_folders{5} = 'Z:\procdata\617m_CPEC\2018-07_IOP2018\dayfiles\fieldversions\vCombined';
sonic_files_folders{6} = 'Z:\procdata\617m_CPEC\2018-08_IOP2018\dayfiles\fieldversions\vCombined';
sonic_files_folders{7} = 'Z:\procdata\617m_CPEC\2018-09_IOP2018\dayfiles\fieldversions\vCombined';
sonic_files_folders{8} = 'Z:\procdata\617m_CPEC\2018-10\dayfiles\fieldversions\vB';
sonic_files_folders{9} = 'Z:\procdata\617m_CPEC\2018-11\dayfiles\fieldversions\vA';
sonic_files_folders{10} = 'Z:\procdata\617m_CPEC\2018-12\dayfiles\fieldversions\vA';
sonic_files_folders{11} = 'Z:\procdata\617m_CPEC\2019-01\dayfiles\fieldversions\vA';
output_folder = 'G:\CPEC_preproc_2019-03-19\tilt_corr';

WINDOW_LENGTH = 18000;	% samples
SONIC_ORIENTATION = 129; % degrees

sector_width = 30;	% width of wind sectors in degrees

mean_uvw = [];
mean_wind_dir = [];
mean_wind_speed = [];
data_intervals_used = 0; % number of data intervals used

% get list of input files
sonic_files = [];
for (n = 1:length(sonic_files_folders))
	sonic_files = [sonic_files; list_files_datestring(sonic_files_folders{n}, {'w'}, '.mat', 'yyyymmdd')];
	fprintf('%s\n', sonic_files_folders{n});
end

fprintf('calc_angle_dependent_tilt_correction: found %d sonic files\n', length(sonic_files));

for (n = 1:length(sonic_files))

	% load sonic data file
	tmp = load([sonic_files(n).path '\\' sonic_files(n).name]);
	%tmpfieldnames = fieldnames(tmp);
	%sonicdata = getfield(tmp, tmpfieldnames{1}); % columns: time, ...
    sonicdata = tmp.data;
	clear tmp tmpfieldnames
	% check data validity and remove invalid data (i.e. where diagnostic word is not zero or wind speed is out of range)
	indices = find(sonicdata(:,7) ~= 0 | sqrt(sonicdata(:,3).^2 + sonicdata(:,4).^2 + sonicdata(:,5).^2) > 30);
	sonicdata(indices, :) = [];

	if (size(sonicdata, 1) > 0)
		for (k = 1:(floor(864000/WINDOW_LENGTH)))
			start_time = round(sonicdata(1,1)) + (k-1)*WINDOW_LENGTH/864000;
			end_time = round(sonicdata(1,1)) + k*WINDOW_LENGTH/864000;
			indices = find(sonicdata(:,1) >= start_time & sonicdata(:,1) < end_time);
			if (length(indices) >= 0.9*WINDOW_LENGTH)
				data_intervals_used = data_intervals_used + 1;
				% calculate mean wind components
				mean_uvw = [mean_uvw; nanmean(sonicdata(indices,3:5), 1)];
				% mean wind speed
				mean_wind_speed = [mean_wind_speed; sqrt(mean_uvw(end,1)^2 + mean_uvw(end,2)^2 + mean_uvw(end,3)^2)];
				% mean wind direction
				mean_wind_dir = [ mean_wind_dir; mod( 180/pi*atan2( -mean_uvw(end,2), mean_uvw(end,1) ) + SONIC_ORIENTATION, 360 ) ];
			end
		end
	end

	progress_bar(length(sonic_files), n);

end

fprintf('%d data intervals processed\n', data_intervals_used);

% sort mean wind data into sectors of wind direction
% 360 overlapping sectors of fixed width
fprintf('sorting data into wind sectors\n');
sorted_data = [];
sorted_data.mean_uvw = [];
sorted_data.mean_wind_speed = [];
sorted_data.mean_wind_dir = [];
sorted_data = repmat(sorted_data, 360, 1);
for (n = 1:length(mean_wind_dir))
	sector_indices = find(abs(mean_wind_dir(n)-(0:359)) < sector_width/2 | abs(mean_wind_dir(n)-360-(0:359)) < sector_width/2 | abs(mean_wind_dir(n)+360-(0:359)) < sector_width/2);
	for (idx = 1:length(sector_indices))
		sorted_data(sector_indices(idx)).mean_uvw = [sorted_data(sector_indices(idx)).mean_uvw; mean_uvw(n,:)];
		sorted_data(sector_indices(idx)).mean_wind_speed = [sorted_data(sector_indices(idx)).mean_wind_speed; mean_wind_speed(n)];
		sorted_data(sector_indices(idx)).mean_wind_dir = [sorted_data(sector_indices(idx)).mean_wind_dir; mean_wind_dir(n)];
	end
	progress_bar(length(mean_wind_dir), n);
end

% calculate tilt correction matrix for each wind sector
fprintf('calculating tilt correction matrices\n');
results.sector_width = sector_width;
results.data_intervals_used = data_intervals_used;
results.data = [];
for (n = 1:360)
	[results.data(n).P, results.data(n).alpha, results.data(n).beta] = calc_tilt_correction(sorted_data(n).mean_uvw(:,1), sorted_data(n).mean_uvw(:,2), sorted_data(n).mean_uvw(:,3), 1);
	results.data(n).count = length(sorted_data(n).mean_wind_dir);
	progress_bar(360, n);
end

% save to file
save([output_folder '\\tilt_correction.mat'], '-struct', 'results');

% plot tilt angles and histogram
alpha = zeros(360, 1);
beta = zeros(360, 1);
count = zeros(360, 1);
for (n = 1:360)
    alpha(n) = results.data(n).alpha;
    beta(n) = results.data(n).beta;
    count(n) = results.data(n).count;
end
subplot1 = subplot(2,1,1);
plot(1:360, alpha*180/pi, 'b', 1:360, beta*180/pi, 'g');
legend('alpha', 'beta', 'Location','southeast');
ylabel('correction [degrees]');
xlim([0 360]);
ylim([-15 15]);
set(subplot1,'XTick',[0 30 60 90 120 150 180 210 240 270 300 330 360]);
subplot2 = subplot(2,1,2);
bar(1:360, count*180/pi);
xlabel('wind direction [degrees]'); ylabel('data count');
xlim([0 360]);
set(subplot2,'XTick',[0 30 60 90 120 150 180 210 240 270 300 330 360]);
