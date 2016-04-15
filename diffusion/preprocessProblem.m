function problemData = preprocessProblem(problemData)
%% Triangulation.
problemData.g = domainSquare(problemData.hmax);
if problemData.isVisGrid,  visualizeGrid(problemData.g);  end
%% Globally constant parameters.
problemData.K = problemData.g.numT;  % number of triangles
problemData.N = nchoosek(problemData.p + 2, problemData.p);  % number of local DOFs
problemData.tau = problemData.tEnd / problemData. numSteps;  % time step size
%% Additional mesh data
problemData.g.markE0Tint = problemData.g.idE0T == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrN = problemData.g.idE0T == 1 | problemData.g.idE0T == 3; % [K x 3] mark local edges on the Neumann boundary
problemData.g.markE0TbdrD = ~(problemData.g.markE0Tint | problemData.g.markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary
problemData.g = computeDerivedGridData(problemData.g);
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
computeBasesOnQuad(problemData.N);
%% Computation of matrices on the reference triangle.
problemData.hatM        = integrateRefElemPhiPhi(problemData.N);
problemData.hatG        = integrateRefElemDphiPhiPhi(problemData.N);
problemData.hatH        = integrateRefElemDphiPhi(problemData.N);
problemData.hatRdiag    = integrateRefEdgePhiIntPhiIntPhiInt(problemData.N);
problemData.hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.N);
problemData.hatSdiag    = integrateRefEdgePhiIntPhiInt(problemData.N);
problemData.hatSoffdiag = integrateRefEdgePhiIntPhiExt(problemData.N);
%% Assembly of time-independent global matrices.
problemData.globM  = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH  = assembleMatElemDphiPhi(problemData.g, problemData.hatH);
problemData.globQ  = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatSdiag, problemData.hatSoffdiag, problemData.g.areaNuE0Tint);
problemData.globQN = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, ...
                      problemData.hatSdiag, problemData.g.areaNuE0TbdrN);
problemData.globS  = problemData.eta * assembleMatEdgePhiPhi(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatSdiag, problemData.hatSoffdiag);
problemData.globSD = problemData.eta * assembleMatEdgePhiIntPhiInt(problemData.g, ...
                      problemData.g.markE0TbdrD, problemData.hatSdiag);
problemData.sysW = [ sparse(2 * problemData.K * problemData.N, 3 * problemData.K * problemData.N) ; ...
                     sparse(problemData.K * problemData.N, 2 * problemData.K * problemData.N), problemData.globM ];
end % function