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

% View the 2D poly in the meshgrid space
polyMatrix = polyObj.matrix2Poly(A);
disp("2D Poly Matrix:");
disp(size(polyMatrix));
disp(polyMatrix);
disp("\n");

% Visualize the 2D poly components
% polyObj.view2DPolyMatrix();
title = "2D Polynomial";
zString = "poly";
polyObj.visAll(polyMatrix, title, zString);

% ---------------------------------------------------------------
% Test the least squares methods for linear and cubic functions
% ---------------------------------------------------------------

% Evaluation pt for the output x vector
xmin = 1; xmax = length(componentNames);
% IDX = randi([xmin, xmax]);  % sample index of sol. vector x
IDX = 1;

% Test a linear function f1
testLinearFunction(polyObj, A, X, Y, IDX);

% Test a cubic function f2
% testCubicFunction(polyObj, A, X, Y, IDX);

% Test the linear function f1(x,y)
function testLinearFunction(polyObj, A, X, Y, IDX)
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
    [kappa, theta, eta] = computeRCNParams(A, b);
    [rcn1, rcn2, rcn3, rcn4] = computeRCNLS(kappa, theta, eta);
    disp("Linear Function F1 RCN Values:");
    disp("RCN1 (y = AA^+b as func of b)     = " + rcn1);
    disp("RCN2 (x^+ = A^+b as func of b)    = " + rcn2);
    disp("RCN3 (y = AA^+b as func of A)     = " + rcn3);
    disp("RCN4 (x^+ = A^+b as func of A)    = " + rcn4);
    disp("=> RCN1 < RCN3 (satisfying Property 3 in NLA)");
    disp("=> RCN4 < RCN2 (satisfying Property 2 in NLA)");
    disp("-------------------------------------------");
    
    % Compute least squares using 3 methods
    computeQRFactorization(polyObj, A, b, IDX);
%     computeNormalEquations(polyObj, A, b, IDX);
%     computeSVD(polyObj, A, b, IDX);
end

% Test the cubic function f2(x,y)
function testCubicFunction(polyObj, A, X, Y, IDX)
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
    [kappa, theta, eta] = computeRCNParams(A, b);
    [rcn1, rcn2, rcn3, rcn4] = computeRCNLS(kappa, theta, eta);
    disp("Cubic Function F2 RCN Values:");
    disp("RCN1 (y = AA^+b as func of b)     = " + rcn1);
    disp("RCN2 (x^+ = A^+b as func of b)    = " + rcn2);
    disp("RCN3 (y = AA^+b as func of A)     = " + rcn3);
    disp("RCN4 (x^+ = A^+b as func of A)    = " + rcn4);
    disp("=> RCN1 < RCN3 (satisfying Property 3 in NLA)");
    disp("=> RCN4 < RCN2 (satisfying Property 2 in NLA)");
    disp("-------------------------------------------");
    
    % Compute least squares using 3 methods
    disp("A:");
    disp(size(A));
    disp(A);
    disp("\n");
    
    computeQRFactorization(polyObj, A, b, IDX);
%     computeNormalEquations(polyObj, A, b, IDX);
%     computeSVD(polyObj, A, b, IDX);
end

% Compute the RCN parameters (kapp, theta, eta)
function [kappa, theta, eta] = computeRCNParams(A, b)
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

% Compute RCN LeastSquare values
function [rcn1, rcn2, rcn3, rcn4] = computeRCNLS(kappa, theta, eta)
    % RCN 1 y as function of b
    rcn1 = 1/cos(theta);
    
    % RCN 2 x as function of b
    rcn2 = kappa/(eta*cos(theta));
    
    % RCN 3 y as function of A 
    rcn3 = kappa/cos(theta);
    
    % RCN 4 x as function of A
    rcn4 = kappa + (kappa^2*tan(theta))/eta;
end

% Compute QR Factorization (Householder Triangulation)
function computeQRFactorization(polyObj, A, b, IDX)
    disp("Householder Triangulation");
    tic
    [Q,R] = qr(A,0);
    x = R\(Q'*b);
    toc
      
    disp("Q:");
    disp("rank(Q) = " + rank(Q));
    disp(size(Q));
    disp(Q);
    disp("\n");
    
    disp("R:");
    disp(size(R));
    disp(R);
    disp("\n");
    
    % Matrix2Components on the Q orthonormal vector set
    QPolyMatrix = polyObj.matrix2Poly(Q);       % Q polynomial matrix
    QRPolyMatrix = polyObj.matrix2Poly(Q*R);    % should equal A
    
    disp("Q Poly Matrix:");
    disp(size(QPolyMatrix));
    disp(QPolyMatrix);
    disp("\n");
    
    % QR Poly Matrix should match the original 2D Poly Matrix
    disp("QR Poly Matrix:");
    disp(size(QRPolyMatrix));
    disp(QRPolyMatrix);
    disp("\n");
    
    % Create the A matrix using the QR decomp
    QR = Q*R;
    err = QR*x - b;
    [M, N] = polyObj.getCoordinates();
    [X, Y] = polyObj.getXYData();
    rErr = reshape(err, M, N);  % reshape the error vector
    disp("Coordinate PPixels:");
    disp("[" + M + ", " + N + "]");
    disp("\n");

    disp("b:");
    disp(b);
    disp("\n");
    disp("x:");
    disp(x);
    disp("\n");
    
    disp("Error = QR*x - b:");
    disp(size(rErr));
    disp(rErr);
    disp("\n");
    disp("X size:");
    disp(size(X));
    disp("Y size:");
    disp(size(X));
    disp("\n");

    % Plot the error in the coordinate system (may want to use visall)
    title1 = "Q Polynomial Matrix";
    title2 = "QR Polynomial Matrix";
    title3 = "Error: Ax - b = (Q*R)*x - b";
    zString1 = "poly";
    zString2 = "error";
    
    polyObj.visAll(QPolyMatrix, title1, zString1);
    polyObj.visAll(QRPolyMatrix, title2, zString1);
    polyObj.visAll(rErr, title3, zString2);
    
    disp(sprintf("QR: x(%d) = %f", IDX, x(IDX)));
    disp("------------------------------------");
    disp("\n");
end

% Compute the normal equations
function computeNormalEquations(polyObj, A, b, IDX)
    disp("Normal Equations:");
    tic
    x = (A'*A)\(A'*b);
    toc
    disp(sprintf("Norm. Eqns: x(%d) = %f", IDX, x(IDX)));
    disp("\n");

    disp("------------------------------------");
end

% Compute the SVD
function computeSVD(polyObj, A, b, IDX)
    disp("SVD:");
    tic
    % Reduced SVD
    [Ur, Sr, Vr] = svd(A, 0);
    r = rank(Sr);
    xr = Vr*(Sr\(Ur'*b));
    Vr = Vr';   % get the original Vr from Vr transpose
    
    % Full SVD
    [Uf, Sf, Vf] = svd(A);
    xf = Vf*(Sf\(Uf'*b));
    Vf = Vf';   % get the original Vf from Vf transpose
    toc
    
    % Reduced SVD
    disp("Sr:");
    disp(size(Sr));
    disp(Sr);
    disp("\n");
    disp("Ur:");
    disp("rank(Ur) = " + rank(Ur));
    disp(size(Ur))
    disp(Ur);
    disp("\n");
    disp("rank(Vr) = " + rank(Vr));
    disp(size(Vr));
    disp(Vr);
    disp("\n");
    
    % Full SVD (1 to r, r+1 to n and equally with rows)
    disp("Sf:");
    disp(size(Sf));
    disp(Sf);
    disp("\n");
    disp("Uf:");
    disp("rank(Uf) = " + rank(Uf));
    disp(size(Uf))
    disp(Uf);
    disp("\n");
    disp("rank(Vf) = " + rank(Vf));
    disp(size(Vf));
    disp(Vf);
    disp("\n");

    
    % Partition Uf into U1 and U2 for first r cols and last m - r
    [m, n] = size(Uf);
    idx1 = 1:r;
    idx2 = (r+1):m;
    U1 = Uf(:,idx1);
    U2 = Uf(:,idx2);
    disp("U1 size:");
    disp(size(U1));
    disp("U2 size:");
    disp(size(U2));
    disp("m - r = " + (m-r));
    disp("U2:");
    disp(U2);
    disp("\n");
    
    % Partition Vf into V1 and V2 first 
    [m, n] = size(Vf);
    if m == r
        disp("Vf FULL RANK!");
        V1 = Vf(:,1:r);
        V2 = V1;
    else
        idx1 = 1:r;
        idx2 = (r+1):m;
        V1 = Vf(:,idx1);
        V2 = Vf(:,idx2);
    end
    disp("rank = " + r);
    disp("V1 size:");
    disp(size(V1));
    disp(V1);
    disp("\n");
    disp("V2 size:");
    disp(size(V2));
    disp(V2);
    disp("\n");
    
    % Create the basis vectors for col. space
    BasisColSpace   = Ur;   % cols of Ur (full rank) form basis of col space
    BasisRowSpace   = Vr;   % cols of Vr (full rank) form basis of row space
    BasisNullSpace  = V2;   % cols of V2 partitioned form basis of null space
    BasisLeftNullSpace = U2;    % cols of U2 form basis of left null space
    
    disp("Basis Col Space:");
    disp("rank = " + rank(BasisColSpace));
    disp(size(BasisColSpace));
    disp(BasisColSpace);
    disp("\n");
    
    disp("Basis Row Space:");
    disp("rank = " + rank(BasisRowSpace));
    disp(size(BasisRowSpace));
    disp(BasisRowSpace);
    disp("\n");
    
    disp("Basis Null Space:");
    disp("rank = " + rank(BasisNullSpace));
    disp(size(BasisNullSpace));
    disp(BasisNullSpace);
    disp("\n");
    
    disp("Basis Left Null Space:");
    disp("rank = " + rank(BasisLeftNullSpace));
    disp(size(BasisLeftNullSpace));
    disp(BasisLeftNullSpace);
    disp("\n");
    
    disp(sprintf("SVD: x(%d) = %f", IDX, xf(IDX)));
    disp("------------------------------------");
    disp("\n");
end
