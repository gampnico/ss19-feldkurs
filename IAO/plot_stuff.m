clear;
clc;
close all;


FFP.f_clim2d = readmatrix("20190613092755310_footprint_raster_fclim2d.txt");
FFP.x_2d = readmatrix("20190613092755310_footprint_raster_x2d.txt");

FFP.y_2d = readmatrix("20190613092755310_footprint_raster_y2d.txt");

surf(FFP(1).x_2d, FFP(1).y_2d, FFP(1).f_clim2d);shading flat;view(2);
hold all;
for i=1:length(r)
    z = FFP(i).fr.*10.*ones(size(FFP(i).yr));
    plot3(FFP(i).xr,FFP(i).yr,z,'r')
end

load