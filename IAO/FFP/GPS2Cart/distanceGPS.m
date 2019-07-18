% Calculates the distance between two GPS points in meters using
% full haversine formula

function d = distanceGPS(latitude1, longitude1, latitude2, longitude2)   
    R = 6371000; % radius of earth in m        
    delta_lat = latitude2-latitude1;
    delta_long = longitude2-longitude1;
    aa = (sind(delta_lat/2))^2 + cosd(latitude1)*cosd(latitude2)*(sind(delta_long/2))^2;
    cc = 2*atan2d(sqrt(aa), sqrt(1-aa));
    d = R*cc*pi/180;
end