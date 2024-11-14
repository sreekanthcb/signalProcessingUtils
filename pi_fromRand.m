%%=================================================================================================
% Synopsis     : Finds value of pi from random numbers
% Last updated : 2021-09-04
%==================================================================================================
clc;clear;close all;

n_points = 1e8; % large value of this gives more accurate value for pi

x       = rand(n_points,1); % array of random numbers between 0 and 1
y       = rand(n_points,1);

z       = sqrt(x.^2+y.^2);% z = (x,y) forms points inside first quadrant of a cartesian plain
cPoints = length(find(z<=1)); % finding how many of this points lies inside a circle of unit radius

CalcPi  = 4*cPoints/n_points; % ratio of area of circle to square is 1/4th of pi (think about it)

fprintf("Calculated value of pi is %f\n",CalcPi);
