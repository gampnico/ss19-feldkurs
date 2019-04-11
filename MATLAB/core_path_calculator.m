%% Initialisation
clear;
clc;
addpath("generated");   % path to output
addpath("topo");    % path to data files

%% Load positional data

% Calculate coordinates from https://epsg.io/31254 
uni = [79638 236722];   % University roof vector
% Receiver vectors:
% hungerburg = [80718 239124]; % Hoehenstrasse 151
% tautermann = [79772 237360]; % Stamser Feld 5
schiesstand = [78881 237443]; % Hotel Schiesstand
receiver_location = schiesstand;

if uni(1)< receiver_location(1)
    x_vector = [uni(1) receiver_location(1)];
else
    x_vector = [receiver_location(1) uni(1)];
end

y_vector = [uni(2) receiver_location(2)];

%% Parse file

% Use preloaded data
load("schiesstand.mat");
name = terrain.name;

% Generate new data
% [xx, yy, zz] = read_tirol_dgm(x_vector, y_vector);
% terrain.name = "Hotel Schiesstand";
% terrain.xx = xx;
% terrain.yy = yy;
% terrain.zz = zz;

%% Get transect
[terrain.path_geo, terrain.path_length, terrain.z_transect] = pathfinder(terrain.xx, terrain.yy, terrain.zz);

%% Plot data

% plotpaths(terrain structure, plot topography, plot transect)
x_geo = plotpaths(terrain, 1, 1);

for i = 1:1:length(x_geo)
    path_function(i) = terrain.z_transect(1) + (terrain.z_transect(end) - terrain.z_transect(1))/(max(x_geo)-min(x_geo)) * x_geo(i);
    path_height(i) = path_function(i) - terrain.z_transect(i);
end
x_geo_normalised = x_geo ./ x_geo(end);
csv_input = [path_height; x_geo_normalised]';

csvwrite("path_height.csv", csv_input)
