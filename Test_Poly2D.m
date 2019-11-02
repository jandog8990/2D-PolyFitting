% ---------------------------------------------
% Test the Implementation of our 2D Poly Class
% ---------------------------------------------

%Poly2D (x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY)
clear all
close all
clc

% Init the low and hi intervals
x_low = -5;
x_hi = 5;
xd = 0.2;   % x step size for x-axis
y_low = -5;
y_hi = 5;
yd = 0.2;   % y step size for y-axis

% 2D poly max degrees
MaxDegreeX = 3;
MaxDegreeY = 3;

% Image dimensionality for 2D poly
M = 10;
N = 5;

% Create a generic 2d poly object
polyObj = Poly2D(x_low, x_hi, xd, y_low, y_hi, yd,... 
    MaxDegreeX, MaxDegreeY, M, N);
% polyObj.components2Matrix();

