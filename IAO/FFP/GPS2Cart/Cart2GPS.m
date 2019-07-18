% This function takes (xx,yy) coordiante (in metres) and the GPS 
% coordinate used as the origin in GPS2Cart.m (in degrees) and converts the
% (xx,yy) coordinate back to GPS, assuming distances are Euclidian 
% (true for distances O(m)) and coordinates not near the poles. The 
% constants 110540 and 111320 are different due to the Earth's oblateness.
% Big H/T to Jim Lewis in this SO thread
% http://stackoverflow.com/questions/2187657/calculate-second-point-knowing-the-starting-point-and-distance

function [lat, long] = Cart2GPS(xx, yy, lat_origin, long_origin)
    delta_long = xx/(111320*cosd(lat_origin));
    delta_lat = yy/110540;

    lat = lat_origin + delta_lat;
    long = long_origin + delta_long;
end