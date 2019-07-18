function [match_index] = overlapCalc(struc1,struc2, overlap)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
match_index = [];
index_array=[];
for n = 1:size(struc2.xn, 1)-1
    for m = 1:size(struc1.xn, 1)-1
        if abs(struc1.xn(1,n) - struc2.xn(1,m)) <overlap
            for i = 1:size(struc2.yn, 2)-1
                for j = 1:size(struc1.yn,2)-1
                    if abs(struc1.yn(i,1) - struc2.yn(j,1))<overlap
                        index_array =[n, m, i, j];
                        match_index = cat(3, match_index, index_array);
                    end
                end
            end
        end
    end
end
end

