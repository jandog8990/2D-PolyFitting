      %% demo() file that demonstrates usage
      function demoSingle()
        close all;
        demoSingleComp();
      end

      function demoSingleComp ()
        % Make a binary image and take a look at the functions!
        N = 18;  % Keep the same dimensions for all of the components.
        M = 18;
        BW        = make_circle(N, M, 6, 5, 4);
        
        % Plot it:
        figure, imagesc(BW), axis image, colormap(gray), title('Single Component Demo');
        
        % Generate coordinate functions from a single binary object        
        imgComp = SingleComp(BW);        
        
        % Standard calls:
        [DistToBdry, BdryLen] = imgComp.CreateBdryCoordSystem();        
        imgComp.PlotBdry ();
        
        % Static calls:
        imgComp.PlotCoords (BW, DistToBdry, BdryLen);                        
      end
 