% ---------------------------------------------
% Test the Implementation of our 2D Poly Class
% ---------------------------------------------

%Poly2D (x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY)
clear all
close all
clc

% Init the low and hi intervals
x_low = -2;
x_hi = 2;
xd = 0.2;   % x step size for x-axis
y_low = -4;
y_hi = 4;
yd = 0.4;   % y step size for y-axis

MaxDegreeX = 3;
MaxDegreeY = 3;

% Create a generic 2d poly object
polyObj = Poly2D(x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY);
polyObj.components2Matrix();

