classdef Components2Matrix < Poly2D
% --------------------------------------------------------
% Class extends the Poly2D class and transfers any 2D 
% poly to the columns of the 2D MatrixForm matrix
% --------------------------------------------------------
  
  	% Need to create a properties file here?
	properties
		poly = []	
		%	n1 = 0
		%	n2 = 0
		%	x = []	% degrees in x direction (row 
		%	y = []
	end

	%% Default methods assume the first arg is an obj
	methods
		
		% Subclass constructor for the Components2Matrix() methods for 
		% transferring the 2D polynomial to columns of 2D MatrixForm matrix
		% 
		% Input:
		%		n1 - degree of the xaxis (row degree)
		%		n2 - degree of the yaxis (col degree)
		%		x - x coords in vector 
		%		y - y coords in vector 
		% TODO: How does this look visually? 1 poly per col? Or for
		% each entry in the Poly matrix we have an associated vector??
		% --------------------------------------------------------------
	  	function compMatrixObj = Components2Matrix(poly, x_low, x_hi, y_low, y_hi, MaxDegreeX, MaxDegreeY)
	  		if (nargin > 0)
				super_args{1} = x_low 
				super_args{2} = x_hi 
				super_args{3} = y_low 
				super_args{4} = y_hi 
				super_args{5} = MaxDegreeX 
				super_args{6} = MaxDegreeY 
			end
			% Call superclass constructor with default args
	  		compMatrixObj@Poly2D(super_args{:})

			% Initialize our object with the poly array
			compMatrixObj.poly = poly
		end
	end
end
