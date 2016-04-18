function [g, NEDGES, NEDNO, NLEDN, NRAEDN, NRIEDN, ETRI, UNRI, UTRI, NSEDN, EMO, EFA] = fort17Read(file, g, NBFR)
%% The program numbers edges of a mesh implicitly, and may have a different edge numbering than specified by the fort.17 file, 
%% but since essential boundary conditions are specified in the fort-files, boundary edges have to be identified and later renumbered.
% evtl TODO: usage in userInput: nur den Teil der Methode verwenden, der
% wirklich gebraucht wird
% TODO: assert g?
assert(isscalar(NBFR) && round(NBFR) == NBFR && NBFR >= 0, 'Invalid number of open sea boundary forcings.');
fileID     = fopen(file);
edgeData   = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
% assert(isequal(round(edgeData), edgeData) && min(edgeData) >= 0, 'Invalid edge information.'); TODO
NEDGES	   = edgeData(1);
assert(NEDGES >= 3, 'Invalid number of edges');
NEDNO      = [edgeData(3:5:5*NEDGES-2), edgeData(4:5:5*NEDGES-1)];
assert(min(NEDNO(:)) >= 1 && max(NEDNO(:)) <= g.numV, 'Non existing nodes appear in edge information.');
NEDEL      = [edgeData(5:5:5*NEDGES  ), edgeData(6:5:5*NEDGES+1)];
assert(min(NEDNO(:)) >= 1 && max(NEDNO(:)) <= g.numT, 'Non existing elements appear in edge information.');
NELED      = [edgeData(1+5*NEDGES+2:4:1+5*NEDGES+4*g.numT-2), ...
				  	  edgeData(1+5*NEDGES+3:4:1+5*NEDGES+4*g.numT-1), ...
					    edgeData(1+5*NEDGES+4:4:1+5*NEDGES+4*g.numT  )];
assert(min(NELED(:)) >= 1 && max(NELED(:)) <= NEDGES, 'Some elements feature non existing edges.');
dataCountr = 1+5*NEDGES+4*g.numT+1;
NIEDS      = edgeData(dataCountr													); dataCountr = dataCountr+1;
NIEDN			 = edgeData(dataCountr+1:2:dataCountr+2*NIEDS -1); dataCountr = dataCountr+2*NIEDS;
if NIEDS  > 0
	assert(min(NIEDN) >=1 && max(NIEDN) <= NEDGES, 'Some interior edge indices exceed the total edge number.');
end % if
NLEDS			 = edgeData(dataCountr													); dataCountr = dataCountr+1;
NLEDN			 = edgeData(dataCountr+1:2:dataCountr+2*NLEDS -1); dataCountr = dataCountr+2*NLEDS;
if NLEDS  > 0
	assert(min(NLEDN) >=1 && max(NLEDN) <= NEDGES, 'Some land edge indices exceed the total edge number.');
end % if
NRAEDS		 = edgeData(dataCountr													); dataCountr = dataCountr+1;
NRAEDN		 = edgeData(dataCountr+1:2:dataCountr+2*NRAEDS-1); dataCountr = dataCountr+2*NRAEDS;
if NRAEDS > 0
	assert(min(NRAEDN) >=1 && max(NRAEDN) <= NEDGES, 'Some radiation edge indices exceed the total edge number.');
end % if
NRIEDS		 = edgeData(dataCountr													); dataCountr = dataCountr+1;
NRIEDN		 = zeros(NRIEDS, 1);
ETRI			 = zeros(NRIEDS, 1);
UNRI			 = zeros(NRIEDS, 1);
UTRI			 = zeros(NRIEDS, 1);
if NRIEDS > 0
	% This procedure is needed, since river edge boundary conditions are epecified in the form 
	% J XI=ETRI(J). U=UNRI(J). V=UTRI(J)
	NRIEDN(1) = edgeData(dataCountr+1);
	edgeData  = fgets(fileID);
	[xi, uv]  = strtok(edgeData);
	[u, v]		= strtok(uv);
	ETRI(1)		= str2double(xi(4:end));
	UNRI(1)		= str2double( u(3:end));
	UTRI(1)		= str2double( v(4:end));
	for J = 2:NRIEDS
		edgeData  = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
		NRIEDN(J) = edgeData(2);
		edgeData	= fgets(fileID);
		[xi, uv]	= strtok(edgeData);
		[u, v]		= strtok(uv);
		ETRI(J)		= str2double(xi(4:end));
		UNRI(J)		= str2double( u(3:end));
		UTRI(J)		= str2double( v(4:end));
	end % for
	assert(min(NRIEDN) >=1 && max(NRIEDN) <= NEDGES, 'Some river edge indices exceed the total edge number.');
	edgeData = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
	assert(isvector(edgeData), 'Invalid open sea boundary edge input.');
	dataCountr = 1;
end % if
NSEDS = edgeData(dataCountr													); dataCountr = dataCountr+1;
NSEDN = edgeData(dataCountr+1:2:dataCountr+2*NSEDS-1); dataCountr = dataCountr+2*NSEDS;
if NSEDS > 0
	assert(min(NSEDN) >= 1 && max(NSEDN) <= NEDGES, 'Some open sea edge indices exceed the total edge number.');
	assert(NBFR >= 1, 'Since there is a open sea boundary, there has to be at least one open sea boundary forcing.');
end % if
%% check for consistency of information
assert( NEDGES == NIEDS + NLEDS + NRAEDS + NRIEDS + NSEDS, ...
				'The total number of edges is not the sum of the number of interior, land, radiation, river and open sea edges.');
assert( isempty(intersect( NIEDN, NLEDN )) && isempty(intersect( NIEDN, NRAEDN)) && ...
				isempty(intersect( NIEDN, NRIEDN)) && isempty(intersect( NIEDN, NSEDN )) && ...
				isempty(intersect( NLEDN, NRAEDN)) && isempty(intersect( NLEDN, NRIEDN)) && ...
				isempty(intersect( NLEDN, NSEDN )) && isempty(intersect(NRAEDN, NRIEDN)) && ...
				isempty(intersect(NRAEDN, NSEDN )) && isempty(intersect(NRIEDN, NSEDN )),   ...
				'Some edges have more than one edge ID.');
EMO = zeros(NBFR, NSEDS);
EFA = zeros(NBFR, NSEDS);
if NSEDS > 0
	for I = 1:NBFR
		EMO(I, :)	 = edgeData(dataCountr  :2:dataCountr+2*NSEDS-2);
		EFA(I, :)	 = edgeData(dataCountr+1:2:dataCountr+2*NSEDS-1);
		dataCountr = dataCountr+2*NSEDS;
	end % for
end % if
fclose(fileID);
end % function
