clear;
clc;
f = readmatrix("q875x.txt");
res = (f(1,end) - f(1,1))/size(f,1);
sim = readmatrix("schiessstand_sim.csv");
q125 =  sim(130, 3);
q250 =  sim(111, 3);
q375 =  sim(93, 3);
q500 =  sim(148/2, 3);
q625 =  sim(56, 3);
q750 =  sim(37, 3);
q875 =  sim(19, 3);
sum_weight = q125 + q250 +q375+q500+q625+q750+q875;

%%
lat = 47.2640351561;
lon = 11.3857068;
load NOx_results_qclass.mat

time_axis = datetime(datevec(results.time));
tl = datetime(2019,05,24,05,00,00);
tu = datetime(2019, 05, 24, 17, 00,00);
tf = isbetween(time_axis, tl, tu);
time_domain = time_axis(tf);
% zm = repmat(results.parameters.SENSOR_HEIGHT, length(time_domain),1);
zm = results.parameters.SENSOR_HEIGHT;
z0 = NaN;
u_mean = results.MET.uvw(tf);
h =1500;
ol = results.MET.L(tf);
windspeed = results.MET.hws(tf);
sigmav = mean(sqrt( (windspeed - u_mean).^2 ),1);
ustar = results.MET.ust(tf);
wind_dir = results.MET.wdir(tf);
wind_dir_mean = mean(wind_dir);
dates = datevec(time_domain);
year = dates(:,1);
month = dates(:,2);
day = dates(:,3);
hours = dates(:,4);
minutes = dates(:,5);

td = length(time_domain);
zm = repmat(zm, length(time_domain),1);
d = repmat(14, length(time_domain),1);
z0 = repmat(-999, length(time_domain),1);
u_mean = repmat(mean(u_mean), length(time_domain),1);
ol = repmat(mean(ol), td, 1);
sigmav = repmat(sigmav, td,1);
ustar = repmat(mean(ustar), td,1);
wind_dir = repmat(0, td,1);
% r = [10:10:80];
% crop = 1;
% csv_input = [
csv_input = [year, day, month, hours, minutes, zm, d, z0, u_mean, ol, sigmav, ustar, wind_dir];
%%
csvwrite("footprint_data.csv", csv_input)