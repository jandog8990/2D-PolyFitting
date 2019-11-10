% ------------------------------------------------
% Poly2D Class
%
% This class takes 2D polynomials and converts
% to matrix form where each column is a poly,
% and vicer versa from matrices to 2D polys
% ------------------------------------------------
classdef Poly2D < handle
    
    properties
        % Rectangular coordinate system variables
        x_low = 0  
        x_hi  = 0
        y_low = 0
        y_hi  = 0
        
        % Pixel sizes along X and Y
        M = 0
        N = 0
        
        % x and y vec for plotting
        xvec = []
        yvec = []
       
        X = []         % X coordinates.
        Y = []         % Y coordinates.
        
        MaxDegreeX = 0 % The maximum polynomial degree for the X-coordinate.
        MaxDegreeY = 0 % The maximum polynomial degree for the Y-coordinate.

        Components = [] % 3D Array. Stores the monomial components:
                         %   1, x, y, x.*y, x.^2, y.^2, etc
                         % in 2D arrays:
                         %   Components(:,:,1) = 1;
                         %   Components(:,:,2) = x; ...
        ComponentNames = ("") % Array os strings that contains the names.
                        % If constructed using "setComponents()" the names
                        % are:
                        %   ComponentNames(1) = "1";
                        %   ComponentNames(2) = "x"; ...
                        % If constructed using "Matrix2Components()" the
                        % names are:
                        %   ComponentNames(1) ="M1";
                        %   ComponentNames(2) ="Mx"; ...
        MatrixForm = [] % Contains the 2D monomials along each column:
                        %  MatrixForm(:,1) = 1; MatrixForm(:,2) = x; ...
    end
    
    %% The default methods assume that the first argument is the object:
    methods              
        %% Main constructor based on the class name.
        % Create the default coordinate system and the 2D polynomials. 

        function imgObj = Poly2D(x_low, x_hi, y_low, y_hi,...
                MaxDegreeX, MaxDegreeY, M, N)
           if (nargin>0)
               setRectCoords(imgObj, x_low, x_hi, y_low, y_hi, M, N);
               setComponents(imgObj, MaxDegreeX, MaxDegreeY);
           end
        end    

        %% Class methods
      
        %% setRectCoords():
        % This function generates a rectangular coordinate system.
        % Inputs:  x_low, x_hi, y_low, y_hi.
        % Outputs: X, Y coordinates stored internally.
        function [] = setRectCoords(imgObj, x_low, x_hi, y_low, y_hi, M, N)
            % Set pixel lengths along axes
            imgObj.M = M;
            imgObj.N = N;
            
            % Setup the low and hi for x vector
            imgObj.x_low = x_low;
            imgObj.x_hi  = x_hi;
            xlen = (abs(x_low) + abs(x_hi) + 1);
            xd = xlen/M;
            imgObj.xvec = x_low:xd:x_hi;
            
            % Setup the low and hi for y vector
            imgObj.y_low = y_low;
            imgObj.y_hi  = y_hi;
            ylen = (abs(y_low) + abs(y_hi) + 1);
            yd = ylen/N;
            imgObj.yvec = y_low:yd:y_hi;
            
            disp("Set Rect Coords:");
            disp("M, N = " + M + ", " + N);
            disp("X vec len = " + length(imgObj.xvec));
            disp("Y vec len = " + length(imgObj.yvec));
            disp("\n");
            
            % Create the mesh grid for poly plots
            [imgObj.X, imgObj.Y] = meshgrid(imgObj.xvec, imgObj.yvec);
        end

        %% setComponents():
        % This function generates a rectangular coordinate system.
        % Inputs:  
        %   MaxDegreeX - x poly max degree
        %   MaxDegreeY - y poly max degree
        %   M - row dimension for image
        %   N - column dimension for image
        % Outputs: Sets up the Components and ComponentNames
        %          as described above.
        function [Z] = setComponents(imgObj, MaxDegreeX, MaxDegreeY)
            % Store the Parameters
            imgObj.MaxDegreeX = MaxDegreeX;
            imgObj.MaxDegreeY = MaxDegreeY;
            X = imgObj.X;
            Y = imgObj.Y;
            disp("Set Components:");
            disp("MaxDegreeX = " + MaxDegreeX);
            disp("MaxDegreeY = " + MaxDegreeY);
            disp("\n");
            
            % Create the component matrices
            syms x y mono
            
            % first layer is the 0th
            b = [];
            n1 = 1;
            n2 = 1;
            disp("Layered 2D Polys:");
            for i = 0:1
                for j = 0:1
                    b = [b x.^i*y.^j];
                    disp(b);
                end
            end
            disp("\n");

            % Loop through the remaining orders and create the matrices
            disp("New B:");
            for n1 = 2:MaxDegreeX
                for n2 = 2:MaxDegreeY
                    % First add the new monomials for x and y for the
                    % current degrees in the loops
                    p1 = x.^n1;
                    p2 = y.^n2;
                    if (~any(b == p1))
                        b = [b p1];
                    end
                    if (~any(b == p2))
                        b = [b p2];
                    end
                    
                    % Loop through remaning orders and create monomials
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
                end
            end
            disp("\n");
            
            % Final check on polynomials
            B = [];
            n1 = MaxDegreeX;
            n2 = MaxDegreeY;
            disp("B:");
            for i = 0:n1
                for j = 0:n2
                    B = [B x.^i*y.^j];
                end
            end
            disp(B);
            disp("\n");
            
            % Create the monomials and put them in vector
            % TODO: Need to keep the order of the poly in the monomials
            % currently we save as y,y^2,x,x^2, etc...
            % HOW???
            p = 0;
            count = 1;
            for i = 0:1:MaxDegreeX
                for j = 0:1:MaxDegreeY
                    p = p + x.^i.*y.^j;
                    mono(count) = x.^i.*y.^j;
                    count = count + 1;
                end
            end
            disp("Final 2D-Poly:");
            disp(p);
            disp("\n");

            disp("Final Monomial Vector:");
            disp(mono);
            disp("\n");
            
            disp("Sorted Monomial:");
            disp(sort(mono));
            disp("\n");
            
            % Create the monomial matrix from the symbolic monomical vector
            % loop through monomials and create 3D matrix
            monoMatrix = [];
            monoNames = [""];
            for i = 1:length(mono)
                m = mono(i);
                monoNames(i) = string(m);
                monoMatrix(:,:,i) = double(subs(m, {x,y}, {X,Y}));
            end

            % Display the contents of the monomial matrix
            disp("Monomial Matrix:");
            for i = 1:1:length(monoNames)
                disp(monoNames(i));
                disp(monoMatrix(:,:,i));
                disp("\n");
            end
            
            % Create the old Z matrix using the symbolic poly
            Zold = double(subs(p, {x,y}, {X, Y}));

            % Sum: loop through all 2D matrices in the 3D mother
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
%             figure();
%             surf(X,Y,Zold);
%             title("Z old using symbols:");
%             
%             figure();
%             surf(X,Y,Znew);
%             title("Z new using matrices:");
        end
        
        % Create 2D Poly by combining 1D X and Y matrices
        function [xyComponents] = create2DPolyComponents(imgObj, N)
        end

        %% components2Matrix
        % Components to matrix method for transferring the 2D polynomial 
        % to columns of 2D MatrixForm matrix
        % 
        % Input:
        %   poly - input polynomial
        %   groupN - the number of polynomial groups
        % TODO: How does this look visually? 1 poly per col? Or for
        % each entry in the Poly matrix we have an associated vector??
        % --------------------------------------------------------------
        function [poly_matrix] = poly2Matrix(imgObj, groupN)           
        end

        %% matrix2Components
        % Matrix to components converts Matrix columns to 2D polys assumes
        % that Coord System and Max Degrees are already defined
        % 
        % Input:
        %   matrix - input 2D matrix with columns of polys
        % TODO: How does this look visually? 1 poly per col? Or for
        % each entry in the Poly matrix we have an associated vector??
        % --------------------------------------------------------------
        function poly = matrix2Poly(imgObj, matrix)
        end
      
    end % End of standard methods
end