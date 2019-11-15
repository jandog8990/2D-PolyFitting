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

        Components = []     % 3D Array. Stores the monomial components:
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
        MatrixForm = []     % Contains the 2D monomials along each column:
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
               setMatrixForm(imgObj);
           end
        end    

        %% Class methods
      
        %% setRectCoords():
        % This function generates a rectangular coordinate system.
        % Inputs:  x_low, x_hi, y_low, y_hi.
        % Outputs: X, Y coordinates stored internally.
        function [] = setRectCoords(imgObj, x_low, x_hi, y_low, y_hi, M, N)
            disp("Set Rect Coords:");
            
            % Set pixel lengths along axes
            imgObj.M = M;
            imgObj.N = N;
            disp("M = " + M);
            disp("N = " + N);
            
            % Setup the low and hi for x vector
            imgObj.x_low = x_low;
            imgObj.x_hi  = x_hi;
            xlen = x_hi - x_low;
            xd = xlen/M;
            x_low = x_low+xd;
            imgObj.xvec = x_low:xd:x_hi;
            disp("Xvec length = " + length(imgObj.xvec));
            
            % Setup the low and hi for y vector
            imgObj.y_low = y_low;
            imgObj.y_hi  = y_hi;
            ylen = y_hi - y_low;
            yd = ylen/N;
            y_low = y_low+yd;
            imgObj.yvec = y_low:yd:y_hi;
            disp("Yvec length = " + length(imgObj.yvec));

            
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
            
            % Create the component matrices
            syms x y mono
            
            % first layer is the 0th
            b = [];
            n1 = 1;
            n2 = 1;
            for i = 0:1
                for j = 0:1
                    b = [b x.^i*y.^j];
                end
            end

            % Loop through the remaining orders and create the matrices
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
                            ex = any(b == p);
                            if (~any(b == p))
                                b = [b p];
                            end
                        end
                    end 
                end
            end
            mono = b;

            disp("Final Monomial Vector:");
            disp("len mono = " + length(mono));
            disp(mono);
            disp("\n");
                        
            % Final check on polynomials
            B = [];
            n1 = MaxDegreeX;
            n2 = MaxDegreeY;
            disp("Algo Monomial Vector:");
            for i = 0:n1
                for j = 0:n2
                    B = [B x.^i*y.^j];
                end
            end
            disp(B);
            disp("\n");
            
            % Create the monomial matrix from the symbolic monomical vector
            % loop through monomials and create 3D matrix
            for i = 1:length(mono)
                m = mono(i);
                ComponentNames(i) = string(m);
                Components(:,:,i) = double(subs(m, {x,y}, {X,Y}));
            end

            % Set the XY Components from the matrix calculations
            imgObj.Components = Components;
            imgObj.ComponentNames = ComponentNames;
        end
        
        % Set the MatrixForm for the 2D polynomial
        function setMatrixForm(imgObj)
            components = imgObj.Components;
            componentNames = imgObj.ComponentNames;
            for i = 1:1:length(componentNames)
                comp = components(:,:,i);
                comp = comp(:);
                imgObj.MatrixForm = [imgObj.MatrixForm comp];
            end
        end
        
        % Visualize 2D Polynomial matrix from components
        function view2DPolyMatrix(imgObj)
            % initialize the system
            syms x y
            X = imgObj.X; Y = imgObj.Y;
            MaxDegreeX = imgObj.MaxDegreeX;
            MaxDegreeY = imgObj.MaxDegreeY;
            Components = imgObj.Components;

            % Sum: loop through all 2D matrices in the 3D mother
            [M,N,P] = size(Components);
            polyMatrix = zeros(M,N);
            for i = 1:1:P
                polyMatrix = polyMatrix + Components(:,:,i);
            end
        end
        
        % Get the X and Y Components matrix
        function Components = getComponents(imgObj)
            Components = imgObj.Components;
        end
        
        % Get the Component names vector (maps to Components matrix)
        function ComponentNames = getComponentNames(imgObj)
            ComponentNames = imgObj.ComponentNames;
        end
        
        % Get the rectangular coords used for the pixel images
        function [M, N] = getCoordinates(imgObj)
            M = imgObj.M;
            N = imgObj.N;
        end
        
        % Get the Vandermonde matrix containing 2D matrix of monomials
        function [vanderMat, componentNames] = getVandermondeMatrix(imgObj)
            % Get the Polynomial Matrix Components (ie. Vandermonde Matrix)
            vanderMat = imgObj.MatrixForm;
            componentNames = imgObj.ComponentNames;
        end
        
        % Get the X and Y input matrices from the meshgrid
        function [X,Y] = getXYData(imgObj)
            X = imgObj.X;
            Y = imgObj.Y;
        end
        
        %% Visualize all components in 2D poly matrix
        function visAll(imgObj, poly, titleString, zString)
            disp("Poly Matrix (size):");
            disp(size(poly));
            disp("X size:");
            disp(size(imgObj.X));
            disp("Y size:");
            disp(size(imgObj.Y));
            disp("\n");
            
            figure()
            surf(imgObj.X, imgObj.Y, poly);
            title(titleString);
            xlabel("X");
            ylabel("Y");
            zlabel("Z (" + zString + ")");
        end
        
        %% components2Matrix
        % Components to matrix method for transferring the 2D polynomial 
        % to columns of 2D MatrixForm matrix
        % 
        % Input:
        %   poly - input polynomial
        %   groupN - the number of polynomial groups
        % --------------------------------------------------------------
        function [poly_matrix] = poly2Matrix(imgObj)           
        end

        %% matrix2Components
        % Matrix to components converts Matrix columns to 2D polys assumes
        % that Coord System and Max Degrees are already defined
        % 
        % Input:
        %   matrix - input 2D matrix with columns of polys
        % --------------------------------------------------------------
        function polyMatrix = matrix2Poly(imgObj, M)
            [m, n] = size(M);
            
            % Loop through the input matrix M and create 3D components
            polyMatrix = zeros(imgObj.M, imgObj.N);
            for i = 1:1:n
                % Reshape the column using given pixel dimensions
                compr = reshape(M(:,i), imgObj.M, imgObj.N);
                
                % Sum the polynomail matrices (i.e. compr sum)
                polyMatrix = polyMatrix + compr;
            end
        end
        function polyMatrix = matrix2Poly2(imgObj, M)
            [m, n] = size(M);
            maxM = imgObj.M^2;
            
            % If the rows don't match the dimension of data append zeroes
            if m ~= maxM
                diffM = maxM - m;
                M = [M; zeros(diffM, n)];
            end
            
            % Loop through the input matrix M and create 3D components
            polyMatrix = zeros(imgObj.M, imgObj.N);
            for i = 1:1:n
                % Reshape the column using given pixel dimensions
                compr = reshape(M(:,i), imgObj.M, imgObj.N);
                
                % Sum the polynomail matrices (i.e. compr sum)
                polyMatrix = polyMatrix + compr;
            end
        end
    end % End of standard methods
end