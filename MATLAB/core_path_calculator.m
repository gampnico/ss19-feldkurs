%% Initialisation
clear;
clc;
addpath("generated");   % path to output
addpath("topo");    % path to data files

%% Load positional data

% Calculate coordinates from https://epsg.io/31254 
uni = [79626 236685];   % University roof vector
% Receiver vectors:
hungerburg = [80722 239124]; % Hoehenstrasse 151
% tautermann = [79772 237360]; % Stamser Feld 5
schiessstand = [78895 237421]; % Hotel Schiesstand
receiver_location = schiessstand;

if uni(1)< receiver_location(1)
    x_vector = [uni(1) receiver_location(1)];
else
    x_vector = [receiver_location(1) uni(1)];
end

y_vector = [uni(2) receiver_location(2)];

%% Parse file

% Use preloaded data
load("hungerburg_dom.mat");
name = terrain.name;

% Generate new data
% [xx, yy, zz] = read_tirol_dgm(x_vector, y_vector, "dom");
% terrain.name = "Gasthaus Schießstand";
% terrain.xx = xx;
% terrain.yy = yy;
% terrain.zz = zz;

%% Get transect
[terrain.path_geo, terrain.path_length, terrain.z_transect] = pathfinder(terrain.xx, terrain.yy, terrain.zz);

%% Plot data

% plotpaths(terrain structure, plot topography, plot transect)
x_geo = plotpaths(terrain, 0, 1);

for i = 1:1:length(x_geo)
    path_function(i) = terrain.z_transect(1) + (terrain.z_transect(end) - terrain.z_transect(1))/(max(x_geo)-min(x_geo)) * x_geo(i);
    path_height(i) = path_function(i) - terrain.z_transect(i);
end
x_geo_normalised = x_geo ./ x_geo(end);
csv_input = [path_height; x_geo_normalised]';
csvwrite("path_height_hungerburg.csv", csv_input)

%% Post Python

sim = readmatrix("hungerburg_sim.csv");

figure(3)
plot(x_geo, path_height, "Color", "k");
hold on
plot(x_geo,sim(:,3)*max(path_height/max(sim(:,3))), "Color", "b")
ylabel("Path height (m)");
xlabel("Distance from UIBK (m)");
title(terrain.name + " Path Height");
legend("Path height", "Path weight", "Location", "northwest")
xlim([0 max(x_geo)])
% text(double(x_geo(end)/2), double(mean(...
%         [max(path_height) min(path_height)],2)), num2str(ceil(mean(path_height))) + "m");