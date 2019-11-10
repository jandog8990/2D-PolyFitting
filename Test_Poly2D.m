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

% Init y value params
y_low = -4;
y_hi = 4;

% 2D poly max degrees
MaxDegreeX = 2;
MaxDegreeY = 2;

% M and N represent pixel lengths along X & Y
M = 5;
N = 5;

% Create a generic 2d poly object
polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
    MaxDegreeX, MaxDegreeY, M, N);

% Test the least squares method
syms x y
c = double(randn());
f = 1 + x + y + c;
disp("Sym function:");
disp(f);
disp("\n");