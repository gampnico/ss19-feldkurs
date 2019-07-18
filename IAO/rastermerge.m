
%% INITIALISATION
clear;
clc;

% Geographic coordinate system
mid.latitude = 47.267123;
mid.longitude = 11.379960;
q3.latitude = 47.268514;
q3.longitude = 11.377360;
q1.latitude = 47.265764;
q1.longitude = 11.382511;

% GPS coordinate system
q3.gps_x = 79029.64;
q3.gps_y = 237155.63;
q1.gps_x = 79423.54;
q1.gps_y = 236855.12;
mid.gps_x = 79228.46;
mid.gps_y = 237003.62;

% read in topographic data
sim = readmatrix("schiessstand_sim.csv");
mid.altitude = sim(148/2, 1);
mid.weight = sim(148/2, 3);
% q3.altitude = 50.1049;
q3.altitude = sim(3*148/4, 1);
q3.weight = sim(3*148/4, 3);
% q1.altitude = 31.4199;
q1.altitude = sim((148 / 4)+1, 1);
q1.weight = sim((148 / 4)+1, 3);

sum_weight = q1.weight + mid.weight + q3.weight;

%% READ IN RASTER DATA

q1.fclim = readmatrix("hq1_raster_fclim2d.txt");
q1.x2d = readmatrix("hq1_raster_x2d.txt");
q1.y2d = readmatrix("hq1_raster_y2d.txt");
q3.fclim = readmatrix("hq3_raster_fclim2d.txt");
q3.x2d = readmatrix("hq3_raster_x2d.txt");
q3.y2d = readmatrix("hq3_raster_y2d.txt");
mid.fclim = readmatrix("mid_raster_fclim2d.txt");
mid.x2d = readmatrix("mid_raster_x2d.txt");
mid.y2d = readmatrix("mid_raster_y2d.txt");

% Fix errors in y2d format
mid.y2dc = repmat(mid.y2d(:, 7), 1, size(mid.x2d, 1));
q1.y2dc = repmat(q1.y2d(:, 7), 1, size(q1.x2d, 1));
q3.y2dc = repmat(q3.y2d(:, 7), 1, size(q3.x2d, 1));

% Weight fclim values and normalise
q1.fclim = q1.weight .* q1.fclim / sum_weight;
mid.fclim = mid.weight .* mid.fclim / sum_weight;
q3.fclim = q3.weight .* q3.fclim / sum_weight;

%% CALCULATE CARTESIAN GRID COORDINATES
% Convert local sub-grid coordinates to underlying Cartesian grid (undergrid),
% based on GPS.
load("raster_data.mat")

% Centrepoint offsets to undergrid axis. Since axes are asymmetric, take larger
% limit to avoid sign errors later.
x_off_q1 = -q1.x2d(1, 1);
x_off_q3 = -q3.x2d(1, 1);
y_off_q1 = -q1.y2dc(1, 1);
y_off_q3 = -q3.y2dc(1, 1);

% Calculate distances between subgrid centrepoints from GPS coordinates.
% Hypotenuse of values should equal ~500m.
x_q3_q1_cd = q1.gps_x - q3.gps_x;
y_q3_q1_cd = q3.gps_y - q1.gps_y;
x_q3_mid_cd = mid.gps_x - q3.gps_x;
y_q1_mid_cd = mid.gps_y - q1.gps_y;

% Calculate centrepoint coordinates relative to undergrid
q1.cp_x = x_off_q3 + x_q3_q1_cd;
mid.cp_x = x_off_q3 + x_q3_mid_cd;
q3.cp_x = x_off_q3;
q1.cp_y = y_off_q1;
q3.cp_y = y_off_q1 + y_q3_q1_cd;
mid.cp_y = y_off_q1 + y_q1_mid_cd;

% Convert local subgrid to undergrid coordinates. Note that q3.yn and mid.yn
% are larger than ylimit.
q3.xn = q3.cp_x + q3.x2d;
q3.yn = q3.cp_y + q3.y2dc;
q1.xn = q1.cp_x + q1.x2d;
q1.yn = q1.cp_y + q1.y2dc;
mid.xn = mid.cp_x + mid.x2d;
mid.yn = mid.cp_y + mid.y2dc;

%% REDIMENSIONALISE
% Grids have differing dimensions. Redimensionalise to smallest grid

q3.xn(end, :) = [];
q3.xn(:, end) = [];
q3.yn(:, end) = [];
q3.yn(end, :) = [];
q3.fclim(:, end) = [];
q3.fclim(end, :) = [];
q3.xn(:, 1) = [];
q3.xn(1, :) = [];
q3.yn(:, 1) = [];
q3.yn(1, :) = [];
q3.fclim(:, 1) = [];
q3.fclim(1, :) = [];

mid.xn(:, end) = [];
mid.xn(end, :) = [];
mid.yn(:, end) = [];
mid.yn(end, :) = [];
mid.fclim(:, end) = [];
mid.fclim(end, :) = [];

% Readjust undergrid origin to account for larger q3
y_off_extra = -q3.yn(1, 1);
q3.yn = y_off_extra + q3.yn;
q1.yn = y_off_extra + q1.yn;
mid.yn = y_off_extra + mid.yn;

o_3_mid = overlapCalc(q3, mid, 1);
o_mid_1 = overlapCalc(mid, q1, 1);
o_3_1 = overlapCalc(q3, q1, 1);

%% SUM OVERLAPS
load("overlap_grid.mat")
ov_3_mid = overlapValues(o_3_mid, q3, mid);
ov_mid_1 = overlapValues(o_mid_1, mid, q1);
% ov_3_1 = overlap_values(o_3_1, q3, q1);

%% GRID MAPPING
load("checkpoint.mat")

% Create underlying fclim grid at 1x1 resolution
% 
% Calculate undergrid limits (Cartesian)
xlimit = x_q3_q1_cd + x_off_q3 + x_off_q1;
ylimit = y_q3_q1_cd + y_off_q1 + y_off_q3 + y_off_extra;
ug.fclim = NaN(ceil(xlimit), ceil(ylimit));

tic
% overwrite undergrid with subgrid fclim values
ug.fclim1 = rasterOverwrite(ug.fclim, q3, 2);
ug.fclim2 = rasterOverwrite(ug.fclim1, mid, 2);
ug.fclim3 = rasterOverwrite(ug.fclim2, q1, 2);
t_elapsed = toc;
%% RASTER OVERWRITE
% overwrite with overlapped values
% ug.fclim4 = ug.fclim3;
% for i = 1:size(ov_3_mid.xy_ol, 1)
%     for j = 1:size(ov_3_mid.xy_ol, 1)
%         
%     end
% end
% ug.fclim5 = rasterOverwrite(ug.fclim4, ov_mid_1, 1);
% t_elapsed_2 = toc;
%% TEST
load("ug_test.mat")
load("cmap_x.mat")
ug.fclim3 = ug.fclim3./sum_weight;
clear Xq Yq Zq ug.fclim4;
ug.fclim4 = ug.fclim3(2000:3500, 1500:3000);
[Xq,Yq] = meshgrid(0:size(ug.fclim4, 2)-1, 0:size(ug.fclim4, 1)-1);
% ug.fclim4(isnan(ug.fclim4))=0;
ug.fclim4(ug.fclim4 < 2e-11) = NaN;
ug.fclim4 = fillmissing(ug.fclim4, "linear");
% Zq = griddata(Xq, Yq, ug.fclim4, Xq, Yq);
% figure(2)
% surf(Xq, Yq, Zq);shading flat;view(2);colormap(cmap_x);
% ug.fclim4 = fillmissing(ug.fclim4, "movmedian", 20);
% ug.fclim5 = fillmissing(ug.fclim4, "nearest");
%%
ug.fclim5 = ug.fclim4;
gaps = find(all(isnan(ug.fclim5),1));
ug.fclim5(:, gaps) = ug.fclim5(:, gaps-1);
ug.fclim5(:, gaps) = ug.fclim5(:, gaps-1);
ug.fclim5(:, gaps) = ug.fclim5(:, gaps-1);
%%
figure(2)
surf(Xq, Yq, ug(1).fclim5);shading flat;view(2);colormap(cmap_x);
colorbar;
caxis([2e-10 5e-6]);
%%
ov_s_mid_1.xn = repmat(ov_mid_1.xy_ol(:,1),1, max(ov_mid_1.idx1(1,1,:)));
ov_s_mid_1.yn = repmat(ov_mid_1.xy_ol(:,1), max(ov_mid_1.idx2(1,1,:)),1);
tic
% ug.fclim6 = rasterOverwrite(ug.fclim4, ov_mid_1, 1);
overlap = 1;
ug.fclim6 = ug.fclim5;
for i = 1:size(ug.fclim5, 1)
    for n =1:length(ov_mid_1.xy_ol)
            if abs(ov_mid_1.xy_ol(n,1) -i)>overlap
                    for j = 1:size(ug.fclim5, 2)
                        if abs(ov_mid_1.xy_ol(n, 2)- j) >overlap
                            ug.fclim6(i,j) = ov_mid_1.mj(n);
                        end
                    end
            end
    end
end
telapse = toc;
%%
tiff = ug.fclim5(:, :).*10^10;
outputFileName = sprintf('smb%d.tiff', 1);
imwrite(tiff,outputFileName,'WriteMode', 'append')
%% INTERPOLATE GRID DATA
% Interpolate subgrid data onto undergrid, overwriting overlaps
% Overwrite overlaps with summed fclim data
%
% % Interpolate q3
%
% [X,Y] = meshgrid(-3:3);
% % equivalent to xn and yn
% V = peaks(X,Y);
% % equivalent to fclim
% [Xq,Yq] = meshgrid(0:round(xlimit), 0:round(ylimit));
% % equivalent to undergrid
%
% Vq = interp2(q3.xn,q3.yn,q3.fclim,Xq,Yq);
% figure(1)
% surf(q3.xn,q3.yn, q3.fclim)
% figure(2)
% surf(Xq,Yq,Vq);
% %%
% % [xq,yq] = meshgrid(0:xlimit, 0:ylimit);
% % vq = griddata(oc_3_mid.xy_ol(:,1), oc_3_mid.xy_ol(:,2), oc_3_mid.mj, xq, yq, "natural");
% % mesh(xq,yq,vq)
% % hold on
% % plot3(oc_3_mid.xy_ol(:,1), oc_3_mid.xy_ol(:,2), oc_3_mid.mj,'o')
% %

% [mi_3_1] = overlap_calc(hq3
% for z = 1:size(match_index, 3)
%     f_ol = hq3raster.fclim()
% end
% for n = 1:length(match_index)
%     mr.overlap_midx(n,:) = midraster.xn(1,n);
%     mr.overlap_midy(:,n) = midraster.yn(n,1);
% end
% %
% % % Expand to subgrid
% % mr.overlap_midx = repmat(mr.overlap_midx, 1, length(match_index));
% % mr.overlap_midy = repmat(mr.overlap_midy, length(match_index), 1);
%
% mr.overlap = cat(3, hq3raster.fclim, midraster.fclim);
%
% for n = 1:length(match_index)
%     idx = match_index(n);
%     mr.fclim(1, idx) = mean(mr.overlap(idx,idx), 3);
% end

%%

% Create raster csv
% V = q3.fclim(:);% csv_input = [];
% V(:,2) = repmat(q3.xn(1,:), 1, length(V)-length(q3.xn(1,:);