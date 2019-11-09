clc
clear
close all

format long
% m is rows n is the polynomial degree
m = 10;
% n = 15;
n = 4;
% t = (0:m-1)'/(m-1);
mv = (0:m-1)';
t = mv/(m-1);
A = [];
for i=1:n
    A = [A t.^(i-1)];
    [M,N] = size(A);
end
disp("Old A:");
disp(A);
disp("\n");

% Temporary vector for testing polynomial order
tmp = A(:,end);
A(:,end) = A(:,2);
A(:,2) = tmp;
disp("New A:");
disp(A);
disp("\n");

b = exp(sin(4*t));
b = b/2006.787453080206;

% disp("b vector:");
% disp(b);
% disp("\n");

x = A\b;
y = A*x;
kappa = cond(A);
theta = asin(norm(b-y)/norm(b));
eta   = norm(A)*norm(b)/norm(y);
disp(sprintf("kappa = %f", kappa));
disp(sprintf("theta = %f", theta));
disp(sprintf("eta   = %f", eta));
disp(" ");

% Compute these after we have figured out the 2D poly funcs
disp("Householder Triangulation")
tic
[Q,R] = qr(A,0);
x = R\(Q'*b);
toc
disp(sprintf("QR: x(%f)= %f", n, x(n)));
disp(" ");

disp("Normal equations");
tic
x = (A'*A)\(A'*b);
toc
disp(sprintf("Norm. Eqns: x(%f)= %f", n, x(n)));

disp("SVD")
tic
[U, S, V] = svd(A, 0);
x = V*(S\(U'*b));
toc
disp(x(n));