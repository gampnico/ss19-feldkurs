function [path_geo, path_length, z_transect] = pathfinder(x, y, z)
%PATHFINDER Returns a linear transect and path parameters between two
% coordinates.
%
% Nicolas Gampierakis 03.2019
%
% Args:
% x: x coordinate vector, MGI coordinates
% y: y coordinate vector, MGI coordinates
% z: height matrix, m
%
% Returns:
% path_geo: geographical path distance, m
% path_length: Euclidean distance, m
% z_transect: height data transect between input coordinates

%% Create transect

[m, n] = size(z);
if m > n, ii = sub2ind([m,n], 1:m, ceil((1:m) * n/m));
elseif m < n, ii = sub2ind([m, n], ceil((1:n) * m/n), 1:n);
else, ii = 1:m+1:m*n;
end

% southerly alignment
z_aligned = flipud(z);
z_transect = z_aligned(ii);
z_transect = z(ii); % if incorrect path
%% Calculate path distances
dx = max(x)-min(x);
dy = max(y)-min(y);
dz = z_transect(end) - z_transect(1);

path_geo = sqrt(dx^2+dy^2);
path_length = sqrt(dz^2+path_geo^2);

end

