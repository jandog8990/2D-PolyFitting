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
xd = 0.1;

% Init y value params
y_low = -4;
y_hi = 4;
yd = 0.2;

% 2D poly max degrees
MaxDegreeX = 2;
MaxDegreeY = 2;

% Create a generic 2d poly object
% polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
%     MaxDegreeX, MaxDegreeY);

syms a b c x
f = a*x^2 + b*x + c;

% C = sym('c', [1 MaxDegreeX*MaxDegreeY+1]);
C = 1:MaxDegreeX*MaxDegreeY+1;
syms x y
p = C(1);
count = 2;
for i = 1:1:MaxDegreeX
    for j = 1:1:MaxDegreeY
        p = p + C(count).*x.^i.*y.^j;
        count = count + 1;
    end
end
disp("Final 2D-Poly:");
disp(p);
disp("\n");

xv = x_low:xd:x_hi;
yv = y_low:yd:y_hi;
disp("len xv = " + length(xv));
disp("len yv = " + length(yv));
            
% Create a simple poly for testing
% [X,Y] = meshgrid(xv, yv);
% Z = subs(p, {x,y}, {X, Y});
% surf(X, Y);

% Create a simple poly for testing
[X,Y] = meshgrid(xv, yv);
Z = double(subs(p, {x,y}, {X, Y}));
% Z = X.*exp(-X.^2 - Y.^2);
surf(X,Y,Z);