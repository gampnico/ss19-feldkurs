function [oc] = overlapValues(overlap_struc, struc1, struc2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for a = 1:size(overlap_struc,3)
    oc.idx1(:,:,a) = overlap_struc(:,1:2:3,a);
    oc.idx2(:,:,a) = overlap_struc(:,2:2:4,a);
    oc.A(a,:) = cat(2,struc1.fclim(oc.idx1(1, 1,a), oc.idx1(1,2,a)), struc2.fclim(oc.idx2(1,1,a), oc.idx2(1,2,a)));
end
oc.mj = sum(oc.A,2);

% xy GPS coordinates 
for a = 1:length(oc.mj)
    oc.xy_ol(a,:) = cat(3,struc1.xn(1, oc.idx1(1, 1,a)), struc1.yn(oc.idx1(1,2,a),1));
end
end

