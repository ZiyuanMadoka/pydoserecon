% MATLAB function
% rotate.m
% Copyright (c) 2019 Ziyuan Wang
% apply rotation matrix with angle (counter-clockwise) in degree to x-y coordinates

function[x_new, y_new] = rotate(x,y,angle) % angle in degree
    x_new = x*cosd(angle) - y*sind(angle);
    y_new = x*sind(angle) + y*cosd(angle);
end