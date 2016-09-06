function problemData = preprocessProblem(problemData)

%% Triangulation.
intBound      = problemData.heightDarcy*ones(problemData.NX+1,1);
problemData.g = generateGridData(width, problemData.NX, problemData.NZdarcy, 1, intBound, intBound+1);

%% Globally constant parameters.
problemData.NT = problemData.g.numTsub;  % number of triangles
problemData.N = (problemData.p + 1)^2;  % number of local DOFs
problemData.tau = problemData.tEnd / problemData. numSteps;  % time step size

%% Additional mesh data
problemData.g.markE0Tint = problemData.g.idE0T(1:problemData.NT,:) == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrN = problemData.g.idE0T(1:problemData.NT,:) == 2 | problemData.g.idE0T(1:problemData.NT,:) == 6; % [K x 3] mark local edges on the Neumann boundary
problemData.g.markE0TbdrD = ~(problemData.g.markE0Tint | problemData.g.markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.NT)

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.p, problemData.ord);

%% Computation of matrices on the reference triangle.
[problemData.hatMc, problemData.hatMx] = integrateRefElemPhiPhi(problemData.p, problemData.ord, problemData.basesOnQuad);
[problemData.hatGc, problemData.hatGx, problemData.hatGy] = integrateRefElemDphiPhiPhi(problemData.p, problemData.ord, problemData.basesOnQuad);
[problemData.hatHc, problemData.hatHx, problemData.hatHy] = integrateRefElemDphiPhi(problemData.p, problemData.ord, problemData.basesOnQuad);
problemData.hatRdiag = integrateRefEdgePhiIntPhiIntPhiInt(problemData.p, problemData.ord, problemData.basesOnQuad);
problemData.hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.p, problemData.ord, problemData.basesOnQuad);
problemData.hatSdiag = integrateRefEdgePhiIntPhiInt(problemData.p, problemData.ord, problemData.basesOnQuad);
problemData.hatSoffdiag = integrateRefEdgePhiIntPhiExt(problemData.p, problemData.ord, problemData.basesOnQuad);

%% Assembly of time-independent global matrices.
problemData.globM  = assembleMatElemPhiPhi(problemData.g, problemData.hatMc, problemData.hatMx);
problemData.globH  = assembleMatElemDphiPhi(problemData.g, problemData.hatHc, problemData.haHx, problemData.hatHy);
problemData.globQ  = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatSdiag, problemData.hatSoffdiag);
problemData.globQN = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, problemData.hatSdiag);
problemData.globS  = problemData.eta * assembleMatEdgePhiPhi(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatSdiag, problemData.hatSoffdiag, problemData.eta);
problemData.globSD = problemData.eta * assembleMatEdgePhiIntPhiInt(problemData.g, ...
                      problemData.g.markE0TbdrD, problemData.hatSdiag, problemData.eta);
problemData.sysW = [ sparse(2 * problemData.NT * problemData.N, 3 * problemData.NT * problemData.N) ; ...
                     sparse(problemData.NT * problemData.N, 2 * problemData.NT * problemData.N), problemData.globM ];
                 
end % function