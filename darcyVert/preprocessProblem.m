function problemData = preprocessProblem(problemData)

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem(1), problemData.numElem(2));

%% Globally constant parameters.
% problemData.NT = problemData.g.numTsub;  % number of triangles
problemData.N = (problemData.p + 1)^2;  % number of local DOFs
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Additional mesh data
problemData.g.markE0Tint = problemData.g.idE0T(1:problemData.NT,:) == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrN = problemData.g.idE0T(1:problemData.NT,:) == 2 | problemData.g.idE0T(1:problemData.NT,:) == 6; % [K x 3] mark local edges on the Neumann boundary
problemData.g.markE0TbdrD = ~(problemData.g.markE0Tint | problemData.g.markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.NT)

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.p, problemData.qOrd);

%% Computation of matrices on the reference triangle.
hatM = integrateRefElemTrapPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatG = integrateRefElemTrapDphiPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatH = integrateRefElemTrapDphiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatRdiag = integrateRefEdgePhiIntPhiIntPhiInt(problemData.p, problemData.qOrd, problemData.basesOnQuad);
hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.p, problemData.qOrd, problemData.basesOnQuad);
hatSdiag = integrateRefEdgePhiIntPhiInt(problemData.p, problemData.qOrd, problemData.basesOnQuad);
hatSoffdiag = integrateRefEdgePhiIntPhiExt(problemData.p, problemData.qOrd, problemData.basesOnQuad);

%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemTrapPhiPhi(problemData.g, hatM);
problemData.globH = assembleMatElemTrapDphiPhi(problemData.g, hatH);
problemData.globQ = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, ...
                      hatSdiag, hatSoffdiag);
problemData.globQN = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, hatSdiag);
problemData.globS = problemData.eta * assembleMatEdgePhiPhi(problemData.g, problemData.g.markE0Tint, ...
                      hatSdiag, hatSoffdiag, problemData.eta);
problemData.globSD = problemData.eta * assembleMatEdgePhiIntPhiInt(problemData.g, ...
                      problemData.g.markE0TbdrD, hatSdiag, problemData.eta);
problemData.sysW = [ sparse(2 * problemData.NT * problemData.N, 3 * problemData.NT * problemData.N) ; ...
                     sparse(problemData.NT * problemData.N, 2 * problemData.NT * problemData.N), problemData.globM ];
                 
end % function