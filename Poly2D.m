classdef Poly2D < handle
    
    properties
        % Rectangular coordinate system variables
        x_low = 0  
        x_hi  = 0
        xd = 0
        y_low = 0
        y_hi  = 0
        yd = 0
        
        xvec = []   % empty data vectors for poly plots
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
        ComponentNames = [] % Array os strings that contains the names.
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

        function imgObj = Poly2D(x_low, x_hi, xd, y_low, y_hi, yd, MaxDegreeX, MaxDegreeY)
           if (nargin>0)
               setRectCoords(imgObj, x_low, x_hi, xd, y_low, y_hi, yd);
               setComponents(imgObj, MaxDegreeX, MaxDegreeY);
           end
        end    

        %% Class methods
      
        %% setRectCoords():
        % This function generates a rectangular coordinate system.
        % Inputs:  x_low, x_hi, y_low, y_hi.
        % Outputs: X, Y coordinates stored internally.
        function [] = setRectCoords(imgObj, x_low, x_hi, xd, y_low, y_hi, yd)
            % Use meshgrid() here.
            
            % Setup the low and hi for x vector
            imgObj.x_low = x_low;
            imgObj.x_hi  = x_hi;
            imgObj.xd = xd;
            
            % Setup the low and hi for y vector
            imgObj.y_low = y_low;
            imgObj.y_hi  = y_hi;
            imgObj.yd = yd;
            
            % Create the xvec and yvec for the mesh grid (plots polys)
            imgObj.xvec = x_low:xd:x_hi;
            imgObj.yvec = y_low:yd:y_hi;
            [imgObj.X, imgObj.Y] = meshgrid(imgObj.xvec, imgObj.yvec);
        end

        %% setComponents():
        % This function generates a rectangular coordinate system.
        % Inputs:  MaxDegreeX, MaxDegreeY
        % Outputs: Sets up the Components and ComponentNames
        %          as described above.
        function [Z] = setComponents(imgObj, MaxDegreeX, MaxDegreeY)        
            % Store the Parameters
            imgObj.MaxDegreeX = MaxDegreeX;
            imgObj.MaxDegreeY = MaxDegreeY;

            % Create a simple poly for testing
            X = imgObj.X;
            Y = imgObj.Y;
            Z = X.*exp(-X.^2 - Y.^2);
            surf(X,Y,Z);
        end

            %% components2Matrix
        % Components to matrix method for transferring the 2D polynomial 
        % to columns of 2D MatrixForm matrix
        % 
        % Input:
        %   poly - input polynomial
        % TODO: How does this look visually? 1 poly per col? Or for
        % each entry in the Poly matrix we have an associated vector??
        % --------------------------------------------------------------
        function [poly_matrix] = components2Matrix(poly)
            disp("Comp 2 mat:")
            disp("poly:")
            disp(poly)
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
        function poly = matrix2Components(matrix)
        end
      
    end % End of standard methods
end