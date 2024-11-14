%%=========================================================================================
% Synopsis     : Finds value of pi from random numbers
%              : idea is ratio of area of circle to square is 1/4th of pi (think about it)
% Last updated : 2021-09-04
%%=========================================================================================
clc;clear;close all;

n_points = 1e8;

x       = rand(n_points,1);
y       = rand(n_points,1);

z       = sqrt(x.^2+y.^2);
cPoints = length(find(z<=1));

CalcPi  = 4*cPoints/n_points;

fprintf("Calculated value of pi is %f\n",CalcPi);
