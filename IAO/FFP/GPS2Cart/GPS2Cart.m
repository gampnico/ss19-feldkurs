% This function takes GPS coordinates and maps them to a Cartesian grid with 
% units of metres. The Cartesian grid is set up so that all the coordiantes
% are in the first quadrant (i.e. all >=0). This is done by first
% defining the origin as the most west and most south points from the GPS
% coordinates, then computing the distance and angle from this origin
% to all the GPS coordinates. The distance is computed using 
% distanceGPS.m  which uses the Haversine formula for a spherical Earth of 
% radius  6371000m. The angle is computed using MATLAB's built-in 
% distance.m function. With the distance and anlge, the x and y coordinates 
% (relative to the defined origin) can be calculated with elementary 
% trigonometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Accuracy is good for coordiantes of the order of a city, but 
% starts to break down for larger distances when the initial 
% bearing != final bearing.


function [x_coord, y_coord, lat_origin, long_origin] = GPS2Cart(lat, long)   

% Identify most west (i.e. min longitude) and most south (i.e. min latitude)
% GPS points (not necesarily the same coordiante) to use as origin for 
% Cartesian system in the first quadrant
% NOTE: can choose different origin but will need to change change
% calculation of x_coord/y_coord accordingly
lat_origin = min(lat);
long_origin = min(long);

% Compute Cartesian coords of GPS data relative to the origin calculated
% above
x_coord = zeros(length(lat) ,1);
y_coord = zeros(length(lat) ,1);
for ii =1:length(lat)    
    hypot = distanceGPS(lat_origin, long_origin, lat(ii), long(ii));
    [~,az] = distance(lat_origin, long_origin, lat(ii), long(ii));
    x_coord(ii) = hypot*sind(az);
    y_coord(ii) = hypot*cosd(az);
end

end