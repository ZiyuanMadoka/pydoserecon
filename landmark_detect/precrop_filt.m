% MATLAB function
% precrop_filt.m
% Copyright (c) 2019 Ziyuan Wang
% find the largest connected component with high intenstiy -- patient
% body, crop the image and apply an n-pixel-by-n-pixel smooth

function [new_origin, I2]=precrop_filt(DRR,npixel)
    
    [last_col,max_comh] = find_comp(DRR,1,0); % column
    [last_row,max_comv] = find_comp(DRR,2,0); % rows
    new_origin = [last_col-max_comh+1, last_row-max_comv+1]; %[column-X, row-Y]
    DRR = DRR(new_origin(2):last_row,new_origin(1):last_col); 
    %imshow(transpose(DRR))
    if npixel >0
        h = ones(npixel,npixel)/npixel^2;  % 3-by-3 pixel averaging filter
        I2 = imfilter(DRR,h);  % average filter
    else
        I2 = DRR;
    end
end
