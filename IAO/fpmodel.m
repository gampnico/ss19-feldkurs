clear;
clc;
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
r = [10:10:80];
crop = 1;
%% Determine domain extent
% Inputs must be scalar
% r = 80
% Preallocate

% x_ci = zeros(1,length(time_domain));
% y_ci = x_ci;
% f_ci = x_ci;
% x_2d = x_ci;
% y_2d = x_ci;
% f_2d = x_ci;
% r = x_ci;
% fr = x_ci;
% xr = x_ci;
% yr = x_ci;

% u_mean_d = mean(u_mean);
% ol_mean = mean(ol);
% ustar_mean = mean(ustar);

% FFP = calc_footprint_FFP(zm, z0, u_mean_d ,h, ol_mean, sigmav, ustar_mean, "r", r, "crop", 1, "rslayer",1, "wind_dir", 0);
% 
td = length(time_domain);
zm = repmat(zm, td, 1);
h = repmat(h, td, 1);
sigmav = repmat(sigmav, td,1);
for n = 1:length(time_domain)
    fieldname = datestr(time_domain(n), "mmmddHHMMSS");
    FFPclim.(fieldname) = calc_footprint_FFP_climatology(zm, z0, u_mean ,h, ol, sigmav, ustar, wind_dir, "r", r, "crop", 1, "rslayer",0, "smooth_data", 1);
    % pass values to structure
end
% 
% fields = fieldnames(FFPclim);
% SX = [FFPclim.(fields{1})];
% 
% for n = 2:length(time_domain)
%     SX = [SX, FFPclim.(fields{n})];
% end
% 
% C = fieldnames(SX);
% SZ = struct(); % scalar structure
% 
% for k = 4:6
%     F = C{k};
%     SZ.(F) = mean(cat(4,SX.(F)),4);
% end
% 
% for n = 1:8
%     SZ(n).r = FFPclim.(fields{1})(n).r;
%     SZ(n).fr = FFPclim.(fields{1})(n).fr;
% end
% % coords are the same for all r
% for k = 6:7
%     for n = 1:8
%         F = C{k};
%         SZ(n).(F) = FFPclim.(fields{1})(n).(F);
%     end
% end
% 
% 
% average values in domain_structure
% 


% j = 1;
% for i= 1:10:60
%   FFP_av(j).field = mean(cat(3, FFP_dom(i:i+9).field), 3); 
%   j = j + 1;
% end
% add values to matrix

% FFP = calc_footprint_FFP_climatology(zm, z0, u_mean, h, ol, sigmav, ustar, wind_dir, "r", r, "crop", 1, "smooth_data", 1);
%h = repmat(1500, length(results.time), 1);
%sigmav = repmat((mean(sqrt( (results.MET.hws - results.MET.uvw).^2 ))),length(results.time),1);
% calculate domain size
%FFP_dom = calc_footprint_FFP(zm,z0,results.MET.uvw,h,results.MET.L,sigmav,results.MET.ust, "r", r, "crop", 1);
% FFP_dom = calc_footprint_FFP(zm,z0,results.MET.uvw,h,results.MET.L,sigmav,results.MET.ust,"wind_dir", results.MET.wdir, "r", r, "crop", 1);


%%
% h = repmat(h, length(time_domain));
% sigmav = repmat(sigmav, length(time_domain));
% 
% FFP = calc_footprint_FFP_climatology(zm, z0, u_mean, h, ol, sigmav, ustar, wind_dir, "r",r, "crop", 1, "smooth_data", 1);
% for n = 1:length(time_domain)
%     fieldname = datestr(time_domain(n), "mmmddHHMMSS");
% FFPclim.(fieldname) = calc_footprint_FFP_climatology(zm, z0, u_mean, h, ol, sigmav, ustar, wind_dir, "r",r, "crop", 1, "smooth_data", 1);
%     % pass values to structure
% end

%load("FFP_day.mat");
fields = fieldnames(FFPclim);
SX = [FFPclim.(fields{1})];
for n = 2:length(time_domain)
    SX = [SX, FFPclim.(fields{n})];
end

C = fieldnames(SX);
SZ = struct(); % scalar structure

for k = 1:3
    F = C{k};
    SZ.(F) = mean(cat(4,SX.(F)),4);
end

for n = 1:8
    SZ(n).r = FFPclim.(fields{1})(n).r;
    SZ(n).fr = FFPclim.(fields{1})(n).fr;
end
% coords are the same for all r
for k = 4:6
    for n = 1:8
        F = C{k};
        SZ(n).(F) = FFPclim.(fields{1})(n).(F);
    end
end
% SZ(n, 8).n = {FFPclim.(fields{1}).n}';

%%
% surf(FFP(1).x_2d, FFP(1).y_2d, FFP(1).fclim_2d);shading flat;view(2);
% hold all;
% for i=1:length(r)
%     z = FFP(i).fr.*10.*ones(size(FFP(i).yr));
%     plot3(FFP(i).xr,FFP(i).yr,z,'r')
% end
%% Crosswind integrated footprint
figure(1)
plot(FFP.x_ci, FFP.f_ci, "k")
hold on
plot([ FFP.x_ci_max FFP.x_ci_max],[0 max(FFP.f_ci)], "b")

%% 2D Plot
figure(2)
surf(SZ(1).x_2d, SZ(1).y_2d, SZ(1).fclim_2d);shading flat;view(2);
hold all;
for i=1:length(r)
    z = SZ(i).fr.*10.*ones(size(SZ(i).yr));
    plot3(SZ(i).xr,SZ(i).yr,z,'r')
end

%% 3D Plot
figure(3)
surf(FFP(1).x_2d, FFP(1).y_2d, FFP(1).f_2d);shading flat