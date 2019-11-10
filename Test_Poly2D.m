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
% polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
%     MaxDegreeX, MaxDegreeY, M, N);

% Test the least squares method
% syms x y
% c = double(randn());
% f = 1 + x + y + c;
% disp("Sym function:");
% disp(f);
% disp("\n");

% Create the x and y test vectors
syms x y;
A = [1];
n1 = 1;
n2 = 3;
disp("A:");
for i = n1:n2
    A = [A x.^i y.^i];
    if (i == 1)
        disp("I = " + i);
        A = [A x.^i*y.^i];
    end
    if (i >= 2)
        disp("I = " + i);
        A = [A x.^i*y.^((i-1):-1:1) x.^((i-1):-1:1)*y.^i x.^i*y.^i];
    end
    disp(A);
end
disp("\n");

% Test polynomial orders
% Should always start with 0 for i,j on the first layer
n1 = 1;
n2 = 1;
b = [];
disp("Layered 2D Polys:");
for i = 0:n1
    for j = 0:n2
        b = [b x.^i*y.^j];
        disp(b);
    end
end
disp("\n");

n1 = 2;
n2 = 2;
b = [b x.^n1 y.^n2];
disp("Layered 2D Polys:");
for i = 1:n1
    for j = 1:n2
        p = x.^i*y.^j;
        disp(p);
        ex = any(b == p);
        disp("exists = " + ex);
        if (~any(b == p))
            b = [b p];
        end
        disp(b);
    end
end
disp("\n");

n1 = 3;
n2 = 2;
p1 = x.^n1;
p2 = y.^n2;
if (~any(b == p1))
    b = [b p1];
end
if (~any(b == p2))
    b = [b p2];
end

disp("Layered 2D Polys:");
for i = 1:n1
    for j = 1:n2
        p = x.^i*y.^j;
        disp(p);
        ex = any(b == p);
        disp("exists = " + ex);
        if (~any(b == p))
            b = [b p];
        end
        disp(b);
    end
end
disp("\n");

B = [];
n1 = 3;
n2 = 2;
disp("B:");
for i = 0:n1
    for j = 0:n2
        B = [B x.^i*y.^j];
    end
end
disp(B);
disp("\n");

