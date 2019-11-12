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
MaxDegreeX = 5;
MaxDegreeY = 2;

% M and N represent pixel lengths along X & Y
M = 100;
N = 5;

% Create a generic 2d poly object
polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
    MaxDegreeX, MaxDegreeY, M, N);

% Get the X and Y components from the matrices
[X, Y] = polyObj.getXYData();
disp("X:");
disp(size(X));
disp("Y:");
disp(size(Y));
disp("\n");

% Get the XY Matrix Components (ie. Vandermonde Matrix)
components = polyObj.getComponents();
componentNames = polyObj.getComponentNames();

% Test the least squares methods

% % Approx. function f1 for testing
% syms x y
% c1 = double(randn());
% f1 = 1 + x + y + c1;
% disp("Function F1:");
% disp(f1);
% disp("\n");
% 
% % Create the f out vector using xvec and yvec
% fout1 = double(subs(f1, {x,y}, {X,Y}));
% b1 = fout1(:);
% disp(sprintf("F1 out size = %s", string(size(fout1))));
% disp(sprintf("b1 out size = %s", string(size(b1))));
% disp("\n");

% Approx. function f2 for testing
syms x y
c2 = double(randn());
f2 = 5 + x^2 + y^3 + c2;
disp("Function F2:");
disp(f2);
disp("\n");

% Create the f out vector using xvec and yvec
fout2 = double(subs(f2, {x,y}, {X,Y}));
b2 = fout2(:);
disp(sprintf("F2 out size = %s", string(size(fout2))));
disp(sprintf("b2 out size = %s", string(size(b2))));
disp("\n");

% Flatten the 2D Component matrices (i.e. 1, x, y, x^2, y^2,...)
% and concatenate in order to form the Vandermonde matrix
vanderMat = [];
for i = 1:1:length(componentNames)
    comp = components(:,:,i);
    comp = comp(:);
    vanderMat = [vanderMat comp];
end
A = vanderMat;
disp("Vandermonde Matrix:");
disp(componentNames);
disp(size(vanderMat));
disp("\n");

% % F1 Approx: solve the system of eqs
% x1 = A\b1;
% y1 = A*x1;
% 
% % Compute RCN, theta and eta
% kappa1 = cond(A);
% theta1 = acos(norm(y1)/norm(b1));
% eta1 = norm(A)*norm(x1)/norm(y1);
% disp(sprintf("kappa1 = %f", kappa1));
% disp(sprintf("theta1 = %f", theta1));
% disp(sprintf("eta1   = %f", eta1));
% disp("\n");

% F1 Approx: solve the system of eqs
x2 = A\b2;
y2 = A*x2;

disp("------------------------------------");

% Compute RCN, theta and eta
kappa2 = cond(A);
theta2 = acos(norm(y2)/norm(b2));
eta2 = norm(A)*norm(x2)/norm(y2);
disp(sprintf("kappa2 = %f", kappa2));
disp(sprintf("theta2 = %f", theta2));
disp(sprintf("eta2   = %f", eta2));
disp("\n");

disp("------------------------------------");

% Evaluation pt for the output x vector
N = 1;

% Compute the QR Orthogonal and Upper Triangular matrices
disp("Householder Triangulation")
tic
[Q,R] = qr(A,0);
x = R\(Q'*b2);
toc

% Show sizes of output vectors Q & R in QR decomp.
% disp("Q Orthog Matrix:");
% disp(size(Q));
% disp("\n");
% 
% disp("R Diag Matrix:");
% disp(size(R));
% disp("\n");

disp("x length = " + length(x));
disp(sprintf("QR: x(%d) = %f", N, x(N)));
disp("\n");

disp("------------------------------------");

% Compute the normal equations
disp("Normal Equations:");
tic
x = (A'*A)\(A'*b2);
toc
disp("x length = " + length(x));
disp(sprintf("Norm. Eqns: x(%d) = %f", N, x(N)));
disp("\n");

disp("------------------------------------");

% Compute the SVD
disp("SVD:");
tic
[U, S, V] = svd(A, 0);
x = V*(S\(U'*b2));
toc

% Show sizes of output vectors USV in SVD decomp.
% disp("U and V sizes:");
% disp(size(U));
% disp(size(V));
% 
% disp("S:");
% disp(size(S));
% % disp(S);
% disp("\n");

disp("x length = " + length(x));
disp(sprintf("SVD: x(%d) = %f", N, x(N)));
disp("\n");
