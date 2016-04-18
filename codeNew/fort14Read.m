function [g, DP, AGRID] = fort14Read(file, ICS, SLAM0, SFEA0)
if nargin == 1
	ICS = 1;
end % if
if nargin <= 2
	SLAM0 = 0; SFEA0 = 0;
end % if
assert(ICS == 1 || ICS == 2, 'Invalid type of coordinate system.');
assert( isscalar(SLAM0) && isscalar(SFEA0), ...
				'Reference coordinates for spherical coordinates (CPP projection) must be real-valued numbers.' );
fileID   = fopen(file  );
AGRID    = fgets(fileID);
gridData = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
fclose(fileID);
NE = gridData(1); 
NP = gridData(2);
assert(round(NE)==NE && NE >= 1, 'Invalid number of elements.');
assert(round(NP)==NP && NP >= 3, 'Invalid number of vertices.');
if ICS == 1
	X		 = gridData(4:4:	4*NP);
	Y		 = gridData(5:4:1+4*NP);
else
	SLAM = gridData(4:4:	4*NP);
	SFEA = gridData(5:4:1+4*NP);
	% Transformation from spherical to cartesian coordinates
	X    = 6378206.4*pi/180*(SLAM-SLAM0)*cos(pi/180*SFEA0);
	Y    = 6378206.4*pi/180*SFEA;
end % if
assert(ismatrix(X) && ismatrix(Y), 'The grid coordinates are invalid.');
DP		  = gridData(6:4:2+4*NP);
assert(ismatrix(DP) && min(DP) >= 0, 'The depth at some vertices is invalid.');
V0Tlist = reshape(gridData(2+4*NP+1:2+4*NP+5*NE), 5, NE).';
assert( isequal(V0Tlist(:,2), 3*ones(NE, 1)), ...
				'Only triangular elements are supported but the input file suggests the usage of a different kind of elements.' );
V0T = V0Tlist(:,3:5);
assert(min(V0T(:)) >= 1 && max(V0T(:)) <= NP, 'Grid connectivity features non existing vertices.');
g.coordV = [X, Y];
g.V0T    = V0T;
g.numT   = NE;
g.numV   = NP;
end % function
