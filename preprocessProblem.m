function problemData = preprocessProblem(problemData)
%% Triangulation.
problemData.g = domainSquare(problemData.hmax);
% Alternative: problemData.g = domainPolygon([0 1 1 0], [0 0 1 1], problemData.hmax);
%% Globally constant parameters.
problemData.K           = problemData.g.numT;  % number of triangles
problemData.N           = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
problemData.tau         = problemData.tEnd / problemData.numSteps;  % time step size

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
%% Computation of matrices on the reference triangle.
problemData.hatM              = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);

problemData.transportData.g = problemData.g;
problemData.transportData.K = problemData.K;
problemData.transportData.tau = problemData.tau;
problemData.transportData.velN = problemData.N;

addpath('transport');
problemData.transportData = preprocessProblem(problemData.transportData);
rmpath('transport');
end % function