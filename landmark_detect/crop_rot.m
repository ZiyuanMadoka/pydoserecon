%MATLAB script
%crop_rot.m
%detect landmarks and retrive measurements from abdominal DRR in mhd format.
% Copyright (c) 2019 Ziyuan Wang
function row_range=crop_rot(Irot)
% find the largest connected component with high intenstiy -- patient
% body and crop the image
    index_zero = Irot>0;
    index_row = sum(index_zero,2);
    row_range = find(index_row>max(index_row)-2); % column
end