% ---------------------------------------------
% Test the Implementation of our 2D Poly Class
% ---------------------------------------------

%Poly2D (x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY)
clear all
close all
clc

% Init the low and hi intervals
x_low = 1;
x_hi = 5;
y_low = 1;
y_hi = 5;

% 2D poly max degrees
MaxDegreeX = 2;
MaxDegreeY = 2;

% Create a generic 2d poly object
polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
    MaxDegreeX, MaxDegreeY);

