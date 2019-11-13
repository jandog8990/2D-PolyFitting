% ---------------------------------------------
% Test the Implementation of our 2D Poly Class
% ---------------------------------------------

%Poly2D (x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY)
clear all
close all
clc

% Init the low and hi intervals
x_low = -1;
x_hi = 1;

% Init y value params
y_low = -2;
y_hi = 2;

% 2D poly max degrees
MaxDegreeX = 2;
MaxDegreeY = 2;

% M and N represent pixel lengths along X & Y
M = 5;
N = 5;

% Create a generic 2d poly object
polyObj = Poly2D(x_low, x_hi, y_low, y_hi, ... 
    MaxDegreeX, MaxDegreeY, M, N);

% Get the X and Y data matrices from meshgrid
[X, Y] = polyObj.getXYData();

% Flatten the 2D Component matrices (i.e. 1, x, y, x^2, y^2,...)
% and concatenate in order to form the Vandermonde matrix
[vanderMat, componentNames] = polyObj.getVandermondeMatrix();
A = vanderMat;
disp("Vandermonde Matrix:");
disp(componentNames);
disp("\n");

% View the 2D poly in the meshgrid space
polyObj.view2DPolyMatrix();

% ---------------------------------------------------------------
% Test the least squares methods for linear and cubic functions
% ---------------------------------------------------------------

% Evaluation pt for the output x vector
N = 1;

% Test a linear function f1
testLinearFunction(A, X, Y, N);

% Test a cubic function f2
testCubicFunction(A, X, Y, N);

% Test the linear function f1(x,y)
function testLinearFunction(A, X, Y, N)
    % Approx. function f1 for testing
    syms x y
    c = double(randn());
    f = 1 + x + y + c;
    disp("Linear Function f:");
    disp(f);
    disp("\n");
    
    % Create the f out vector using xvec and yvec
    fout = double(subs(f, {x,y}, {X,Y}));
    b = fout(:);
    disp("------------------------------------");
    
    % Compute RCN params for RCN values
    computeRCNParams(A, b);
    
    % Compute least squares using 3 methods
    computeQRFactorization(A, b, N);
    computeNormalEquations(A, b, N);
    computeSVD(A, b, N);
end

% Test the cubic function f2(x,y)
function testCubicFunction(A, X, Y, N)
    % Approx. function f2 for testing
    syms x y
    c = double(randn());
    f = 5 + x^2 + y^3 + c;
    disp("Cubic Function f:");
    disp(f);
    disp("\n");

    % Create the f out vector using xvec and yvec
    fout = double(subs(f, {x,y}, {X,Y}));
    b = fout(:);
    disp("------------------------------------");
    
    % Compute RCN params for RCN values
    computeRCNParams(A, b);
    
    % Compute least squares using 3 methods
    computeQRFactorization(A, b, N);
    computeNormalEquations(A, b, N);
    computeSVD(A, b, N);
end

% Compute the RCN parameters (kapp, theta, eta)
function computeRCNParams(A, b)
    % f Approx: solve the system of eqs
    x = A\b;
    y = A*x;

    % Compute RCN, theta and eta
    kappa = cond(A);
    theta = acos(norm(y)/norm(b));
    eta = norm(A)*norm(x)/norm(y);
    disp(sprintf("kappa = %f", kappa));
    disp(sprintf("theta = %f", theta));
    disp(sprintf("eta   = %f", eta));
    disp("------------------------------------");
    disp("\n");
end

% Compute QR Factorization (Householder Triangulation)
function computeQRFactorization(A, b, N)
    disp("Householder Triangulation");
    tic
    [Q,R] = qr(A,0);
    x = R\(Q'*b);
    toc
    disp(sprintf("QR: x(%d) = %f", N, x(N)));
    disp("------------------------------------");
    disp("\n");
end

% Compute the normal equations
function computeNormalEquations(A, b, N)
    disp("Normal Equations:");
    tic
    x = (A'*A)\(A'*b);
    toc
    disp(sprintf("Norm. Eqns: x(%d) = %f", N, x(N)));
    disp("\n");

    disp("------------------------------------");
end

% Compute the SVD
function computeSVD(A, b, N)
    disp("SVD:");
    tic
    [U, S, V] = svd(A, 0);
    x = V*(S\(U'*b));
    toc
    disp(sprintf("SVD: x(%d) = %f", N, x(N)));
    disp("------------------------------------");
    disp("\n");
end
