function [output_fclim] = rasterOverwrite(undergrid_fclim, raster, overlap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
output_fclim = undergrid_fclim;
for n = 1:size(undergrid_fclim, 1) - 1
    for m = 1:size(raster.xn, 1) - 1
        if abs(raster.xn(1, m) - n) < overlap
            for i = 1:size(undergrid_fclim, 2)-1
                for j = 1:size(raster.yn,2)-1
                    if abs(raster.yn(j,1) - i) < overlap
                        output_fclim(n, i) = raster.fclim(m, j);
%                     else
%                         output_fclim(n, i) = undergrid.fclim(n, i);
                    end
                end
            end
        end
    end
end
end

