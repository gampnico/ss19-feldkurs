function [x_geo] = plotpaths(terrain, topo_arg, trans_arg)
%PLOTPATHS Plots topography and path transect.
% Nicolas Gampierakis 03.2019
%
% Args:
% terrain: terrain structure variable.
% topo_arg: switch value between 0 and 1 to plot topographical map.
%trans_arg: switch value between 0 and 1 to plot transect.

% Returns:
% topo_map: surface plot of x,y,z.
% transect_plot: plot of path transect.

%% Topographical map
if topo_arg == 1
    figure(1)
    topo_map = surf(terrain.xx, terrain.yy, terrain.zz);
    % topo_map.EdgeColor = [47/255 79/255 79/255];
    topo_map.EdgeColor = "none";
    demcmap(terrain.zz)
    topo_map.FaceColor = "texturemap";
    topo_map.FaceLighting = "flat";

    % Prettify
    % Divide the lengths by the number of lines needed
    xnumlines = 0; % 10 lines
    ynumlines = 110; % 10 partitions
    xspacing = round(length(terrain.xx) / xnumlines);
    yspacing = round(length(terrain.yy) / ynumlines);
    % Plotting lines in the X-Z plane
    hold on
    for i = 1:yspacing:length(terrain.yy)
        Y1 = terrain.yy(i) * ones(size(terrain.xx)); % a constant vector
        Z1 = terrain.zz(i, :);
        plot3(terrain.xx, Y1, Z1, '-k');
    end
    % Plotting lines in the Y-Z plane
    for i = 1:xspacing:length(terrain.xx)
        X2 = terrain.xx(i) * ones(size(terrain.yy)); % a constant vector
        Z2 = terrain.zz(:, i);
        plot3(X2, terrain.yy, Z2, '-k');
    end
    hold off
end

%% Transect map
if trans_arg == 1
    figure(2);
    terrain.z_transect=fliplr(terrain.z_transect);
    %normalise path distance when plotting
    x_geo = 1:terrain.path_geo / length(terrain.z_transect):terrain.path_geo;
    plot(x_geo, terrain.z_transect, "Color", "k");
    line([x_geo(1) x_geo(end)], [terrain.z_transect(1) terrain.z_transect(end)],...
        "Color", "b", "LineStyle", "--");
    ylabel("Altitude (m)");
    xlabel("Distance from UIBK (m)");
    title(terrain.name);
    text(double(x_geo(end)/2), double(mean(...
        [max(terrain.z_transect) min(terrain.z_transect)],2)), num2str(ceil(terrain.path_length)) + "m");
    xlim([0 max(x_geo)])
end

end

