% ---------------------------------------------
% Test the Implementation of our 2D Poly Class
% ---------------------------------------------

%Poly2D (x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY)
clear all
close all
clc

% Init the low and hi intervals
x_low = -2;
xd = 1;
x_hi =2;

% Init y value params
y_low = -4;
yd = 1.5;
y_hi = 4-yd;

% 2D poly max degrees
MaxDegreeX = 1;
MaxDegreeY = 1;

% Create a generic 2d poly object
% polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
%     MaxDegreeX, MaxDegreeY);

syms a b c x
f = a*x^2 + b*x + c;

% C = sym('c', [1 MaxDegreeX*MaxDegreeY+1]);
% C = 1:(MaxDegreeX+1)*(MaxDegreeY+1);
C = ones(1, (MaxDegreeX+1)*(MaxDegreeY+1));
syms x y mono
p = 0;
count = 1;
for i = 0:1:MaxDegreeX
    for j = 0:1:MaxDegreeY
        p = p + C(count).*x.^i.*y.^j;
        mono(count) = x.^i.*y.^j;
        disp(p);
        count = count + 1;
    end
end
disp("Final 2D-Poly:");
disp(p);
disp("\n");

disp("Final Monomial Vector:");
disp(mono);
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
[M,N] = size(X);    % need the size of the X input to represent X
Zold = double(subs(p, {x,y}, {X, Y}));

% Create the monomial matrix from the symbolic monomical vector
% monoMatrix = double(subs(mono, {x,y}, {X,Y}));
% loop through monomials and create 3D matrix
monoMatrix = [];
monoNames = [""];
for i = 1:length(mono)
    m = mono(i);
    monoNames(i) = string(m);
    monoMatrix(:,:,i) = double(subs(m, {x,y}, {X,Y}));
end

disp("Monomial Matrix:");
for i = 1:1:length(monoNames)
    disp(monoNames(i));
    disp(monoMatrix(:,:,i));
    disp("\n");
end

[M,N,P] = size(monoMatrix);
Znew = zeros(M,N);
for i = 1:1:P
    Znew = Znew + monoMatrix(:,:,i);
end

disp("Z old:");
disp(Zold);
disp("\n");

disp("Z New:");
disp(Znew);
disp("\n");

disp("Matrix Comparison (Zdiff):");
Zdiff = Znew - Zold;
disp(Zdiff);
disp("\n");

% Display the final Z matrices for the old way and new way
figure();
surf(X,Y,Zold);
title("Z old using symbols:");

figure();
surf(X,Y,Znew);
title("Z new using matrices:");
