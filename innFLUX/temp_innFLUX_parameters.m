% path definitions
output_folder = 'G:\IAO2018_PTR\IAO2018_JAS';
sonic_files_folders{1} = 'Z:\procdata\617m_CPEC\2018-07_IOP2018\dayfiles\fieldversions\vCombined';
sonic_files_folders{2} = 'Z:\procdata\617m_CPEC\2018-08_IOP2018\dayfiles\fieldversions\vCombined';
sonic_files_folders{3} = 'Z:\procdata\617m_CPEC\2018-09_IOP2018\dayfiles\fieldversions\vCombined';

tracer_files_folders = {};  %{1} = 'G:\IAO2018_PTR\IAO-preIOP2018\NO_mode\vC\concentrations';
%tracer_files_folders{2} = 'G:\IAO2018_PTR\IAO-IOP2018_Summer\NO_mode\vC\concentrations';
%tracer_files_prefix = 'ptrms';
%tracer_files_folders{1} = 'C:/FluxDATA/CLD/2015';
%tracer_files_prefix = 'cld';
tracer_files = {};
%tracer_files{1} = 'C:/FluxDATA/PTRMS/PTRMS_toluene.mat';
%tracer_files{1} = 'C:/FluxDATA/CLD/2015/cld20150701_20151020_5Hz.mat';
%tracer_files{1} = 'C:/FluxDATA/CLD/2016/CLD_20160513-45.mat';
%tracer_files{1} = 'C:/FluxDATA/CaRDS/CaRDS_20160513-45.mat';
tilt_correction_filepath = 'Z:\procdata\617m_CPEC\IOP2018_Summer\tilt_corr\tilt_correction.mat';
pressure_filepath = 'Z:\procdata\617m_CPEC\IOP2018_Summer\dayfiles\vH\aux_and_1min\pressure.mat';

% IRGA data: if sonic files contain IRGA data, define column indices here, otherwise set them to []
irga_columns = [8 9];
irga_names = {'CO2', 'H2O'};
irga_H2O_index = 2; % needed for WPL correction; determines which of the IRGA tracers is H2O
irga_flag_column = 10;

% tracers to process (use [] to process ALL tracers)
tracer_indices = [];
%tracer_indices = [94 96 126 127 204 235];
%tracer_indices = [94 96 126 127];
%tracer_indices = [56 126];

% data segments, for which a global lagtime is used each
%data_segments = [datenum(2015,7,1) datenum(2015,7,21) datenum(2015,8,7) datenum(2015,10,21)];
data_segments = [];

% parameters
SONIC_ORIENTATION = 129;        % degrees
SENSOR_HEIGHT = 24.86;          % sensor height in meters above roughness height
WINDOW_LENGTH = 18000;          % samples
SAMPLING_RATE_SONIC = 10;       % per second; sonic anemometer sampling rate
SAMPLING_RATE_TRACER = 10;      % per second; tracer concentration sampling rate
DISJUNCT_EC = 0;                % if 1, disjunct eddy covariance is applied
DETREND_TRACER_SIGNAL = 1;      % if 1, tracer signal detrending is applied
MAX_LAG = 2000;                 % samples; maximum lag for covariance function
LAG_SEARCH_RANGE = 100;         % samples; maximum for lag-time search (+/-)
COVPEAK_FILTER_LENGTH = 15;     % samples; smoothing length of filter applied prior to finding covariance peak
NUM_FREQ_BINS = 60;             % number of logarithmically spaced frequency bins for cospectra
FREQ_BIN_MIN = 1e-3;            % lowest frequency bin for scaled cospectra
FREQ_BIN_MAX = 1000;            % highest frequency bin for scaled cospectra
APPLY_TILT_CORRECTION = 1;      % if 1, sonic tilt correction is applied
APPLY_WPL_CORRECTION = 0;       % if 1, WPL correction is applied
COMPLETENESS_THRESHOLD = 0.5;   % threshold value for completeness of input data (0.0 - 1.0)
DEFAULT_PRESSURE = 1000;        % hPa / mbar; fixed pressure used if no pressure is given in the input files
