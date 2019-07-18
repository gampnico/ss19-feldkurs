% Demonstration of GPS2Cart.m function which takes GPS coordinates and
% maps them to the first quadrant of a Cartesian grid (in units of metres).
% Also demonstrates Cart2GPS.m function which takes takes a GPS origin and 
% maps Cartesian coordiantes back to GPS.
%
% By Tom Ashbee
% Requires  GPS2Cart.m, Cart2GPS.m, distanceGPS.m, plot_google_map.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Accuracy is good for coordiantes of the order of a city, but 
% starts to break down for larger distances when distances are not
% Euclidian.


clear

%% GPS data 
% Emirates Stadium, White Heart Lane, Stamford Bridge, Upton Park
lat = [51.5550 51.6033 51.4817 51.5319];
long = [-0.1086 -0.0658 -0.1911 0.0394];

%% Plot GPS data and a Google map
figure(4)
clf
subplot(2,1,1)
plot(long,lat,'.r','MarkerSize',20);
hold on
plot_google_map
ylabel({'$\phi$ [degrees]'},'interpreter','latex','FontSize',20)
xlabel({'$\lambda$ [degrees]'},'interpreter','latex','FontSize',20)


%% Convert GPS data to Cartesian data
[x_coord, y_coord, lat_origin, long_origin] = GPS2Cart(lat, long);

%% Plot Cartesian data
subplot(2,1,2)
plot(x_coord, y_coord,'.b','MarkerSize',20)
ylabel({'$y$ [metres]'},'interpreter','latex','FontSize',20)
xlabel({'$x$ [metres]'},'interpreter','latex','FontSize',20)
axis equal
axis([0 max(x_coord) 0 max(y_coord)])


%% Convert Cartesian data back to GPS data
% If nothing else, this is useful to test how good approximations are
[lat, long] = Cart2GPS(x_coord, y_coord, lat_origin, long_origin)


