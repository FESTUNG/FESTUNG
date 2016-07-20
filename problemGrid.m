function [g, markE0TL, markE0Taux, depth] = problemGrid(problem, structure, initialSize, refinement, isVisGrid)
%% define the grid and the boundary edges where land boundary conditions should hold
switch problem
  case 1
    X = [0 100 100 0]; Y = [0 0 100 100];
    switch structure
      case 0 % unstructured grid
        g = domainHierarchy(X, Y, initialSize, refinement);
      case 1 % Friedrichs-Keller triangulation
        g = domainSquare(0.5^refinement);
      otherwise
        error('Unknown grid option.');
    end % switch
    %% boundary specification, the user must specify where the land boundary should be, the rest of the edge identification is done by the program itself
    markE0TL  = g.idE0T == 1 | g.idE0T == 3 ... 
                             | g.idE0T == 4;
    depth = 0; % not needed
  case 2 % bahamas
%     [coordV , V0T,   depth] = readGridData  ();
%     [OSNodes, LNodes      ] = readBoundaries();
%     g = domainBuild(coordV, V0T, OSNodes, LNodes);
%     markE0TL = g.idE0T == 1 | g.idE0T == 2;
		[g depth] = fort2Mat('bahamas');
		 markE0TL = g.idE0T == 1;
  otherwise
    error('Unknown grid.');
end
markE0Taux = cell(3,3); % auxiliary cell of vectors of length K, needed in some routines
for nn = 1:3
  for np = 1:3
    markE0Taux{nn,np} = g.markE0TE0T{nn,np} * ones(g.numT,1);
  end
end
if isVisGrid
  visualizeGrid(g);
end % if
end % function