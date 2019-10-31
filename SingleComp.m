classdef SingleComp < handle
    
    properties
        Comp = []  % Black and white image
        
        X = []         % X boundary coordinates.
        Y = []         % Y boundary coordinates.
        
        Xcen = 0    % X centroid
        Ycen = 0    % Y centroid
    end
    
    %% The default methods assume that the first argument is the object:
    methods              
       %% Main constructor based on the class name.
       % Create the boundary from a component. 
       function imgObj = SingleComp (Comp)
           if (nargin>0)
               % Data starts at the second entry:
               C = contourc (double(Comp), [0.5, 0.5]);
               imgObj.X = C(1, 2:end);
               imgObj.Y = C(2, 2:end);  
               imgObj.Comp = Comp;
               
               % Compute the centroid:
               s  = regionprops(Comp, 'centroid');
               imgObj.Xcen = s.Centroid(1);
               imgObj.Ycen = s.Centroid(2);   
           end
       end    
      
      %% Distance functions for a single component:
      
        %% CreateBdryCoords():
        %   Creates the boundary coordinate system from X, Y.        
        % Outputs: (DistToBdry, BdryLen): The coordinates of the distances to the boundary.
      function [DistToBdry, BdryLen] = CreateBdryCoordSystem (imgObj)        
        BdryLength = BdryToLength (imgObj);                      % Get the length of the boundary.
        [DistToBdry, BdryInd] = DistCompToBdry (imgObj);  % Get the distance to the boundary.
        BdryLen = BdryLength(BdryInd);
      end
      
        %% BdryToLength():
        %   Computes the length of the CLOSED boundary from the first point to the length of the array.        
        % Outputs: BdryLength: holds the lengths along the array.
      function [BdryLength] = BdryToLength (imgObj)        
        % Take distances as you go:
        len = length(imgObj.X);
        DistPoints  = sqrt( (imgObj.X(2:len) - imgObj.X(1:len-1)).^2 + ...
                                     (imgObj.Y(2:len) - imgObj.Y(1:len-1)).^2 );
        BdryLength = [0.0, cumsum(DistPoints)];
      end

        %% DistCompToBdry():
        %     Computes all of the distances from each point to the boundary of a binary component.        
        % Outputs:
        %     DistToBdry:    Signed distance to the boundary:      
        %     BdryInd:         (BdryX [BdryInd], BdryY [BdryInd]) are the coordinate points that lie on the
        %                           boundary and they are the closest to each point in the boundary.
      function [DistToBdry, BdryInd] = DistCompToBdry (imgObj)
        % Generate a simple coordinate system:
        N = size(imgObj.Comp, 1);     % Number of rows -> y
        M = size(imgObj.Comp, 2);     % Number of columns -> x

        % Generate an image of the same dimensions
        DistToBdry = zeros(N, M);
        BdryInd      = zeros(N, M);

        % Compute the distances to the boundary from each point in the image:
       len = length(imgObj.X);
       for x = 1:M
         for y=1:N
           Dist = sqrt((imgObj.X - x).^2 + (imgObj.Y - y).^2 );
           [minDist, Ind] = min (Dist);        
           DistToBdry (x, y) = minDist;
           BdryInd      (x, y) = Ind;
         end
       end

      end

        %% PlotBdry() shows the boundary of the component
       function PlotBdry (imgObj)
          % Plot them
          figure;
          imagesc(imgObj.Comp), axis image, colormap(gray); hold on;
          contour(imgObj.Comp, [0.5, 0.5]);

          % Plot them again
          figure;
          imagesc(imgObj.Comp), axis image, colormap(gray); hold on;
          plot(imgObj.X, imgObj.Y, '*');
          plot([imgObj.X, imgObj.X(1)], [imgObj.Y, imgObj.Y(1)]);
       end

    end % End of standard methods
    
   %% Static methods do not require anything!
   % We are simply using the namespace.
    methods (Static)               
       %% Plot functions
       function PlotCoords (Comp, DistToBdry, BdryLen)
          figure
          subplot(1,3,1), imagesc(Comp), axis image, colormap(gray), title('Binary Component');
          subplot(1,3,2), imagesc(DistToBdry), axis image, colormap(gray), title('Dist. to Boundary');
          subplot(1,3,3), imagesc(BdryLen), axis image, colormap(gray), title('Bdry Length');
       end       
       
    end % Static Methods
    
end