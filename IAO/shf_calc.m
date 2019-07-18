clear;
clc;
load NOx_results_qclass.mat

r_d = 287.05;
cp = 1005.7;

rho_air = results.MET.p*100 ./ (results.MET.T .* r_d);
shf = results.MET.wtheta .* cp .* rho_air;
ws = results.MET.uvw;
ust = results.MET.ust;

lob= results.MET.L;
time = datetime(datevec(results.time));
tl = datetime(2019,05,24,03,00,00);
tu = datetime(2019, 05, 24, 15, 56,00);
tf = isbetween(time, tl, tu);
t_axis = time(tf);
shf = shf(tf);
ws = ws(tf);
lob = lob(tf);
ust = ust(tf);

opts = detectImportOptions('ffp_export2.csv');
ffp = readmatrix("ffp_export2.csv", opts);
shf_exp = ffp(:,5);

date_exp = ([ffp(:,8), ffp(:,9), ffp(:,10), ffp(:,11), ffp(:,12), zeros(length(ffp), 1)]);
date_exp = datetime(date_exp);
% figure(1)
% plot(t_axis, shf)
% hold on
% plot(date_exp, shf_exp)
% 
n = 26;
s1 = size(shf_exp, 1);      % Find the next smaller multiple of n
m  = s1 - mod(s1, n);
y  = reshape(shf_exp(1:m), n, []);     % Reshape x to a [n, m/n] matrix
shf_30 = transpose(sum(y, 1) / n);  % Calculate the mean over the 1st dim

% figure(2)
% plot(t_axis,shf_30)
% hold on
% plot(t_axis, shf)
% scatter(shf_30, shf)
% shf_int = interp(shf, 27);
% shf_int(end-1) = [];
figure(2)
plot(t_axis, shf)
% figure(3)plot(t_axis, shf)

% plot(date_exp, shf_exp)
% hold on
% plot(date_exp, shf_int)
% figure(4)
% scatter(shf_int, shf_exp)
dv = datevec(t_axis);
%%
plot(t_axis,ws)
csv_input = [dv(:, 1), dv(:,2), dv(:,3), dv(:, 4), dv(:, 5), dv(:, 6), shf, ws, lob];
csvwrite("rooftop_hungerburg.csv", csv_input)