function problemData = preprocessProblem(problemData)

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Additional mesh data
% [K x 3] mark local edges that are interior or boundary
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g); 
problemData.g.markE0TbdrN = sparse(problemData.generateMarkE0TbdrN(problemData.g));
problemData.g.markE0TbdrD = sparse(problemData.generateMarkE0TbdrD(problemData.g));
problemData.g.markE0TbdrCoupling = sparse(problemData.generateMarkE0TbdrCoupling(problemData.g));

%% Configuration output.
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('Running testcase "%s".\n', problemData.testcase);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d x %d (%d) trapezoids.\n', ...
        problemData.p, problemData.N, problemData.numElem(1), problemData.numElem(2), problemData.g.numT);
if problemData.isStationary
  fprintf('Computing stationary solution.\n');
else
  fprintf('%d time steps from t = %g to %g.\n', problemData.numSteps, problemData.t0, problemData.tEnd);
end % if
fprintf('-------------------------------------------------------------------------------------------\n');

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuadTensorProduct(problemData.p, struct, [problemData.qOrd, problemData.qOrd+1]);

%% Computation of matrices on the reference element.
problemData.hatM = integrateRefElemTetraPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatG = integrateRefElemTetraDphiPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
hatH = integrateRefElemTetraDphiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatRdiag = integrateRefEdgeTetraPhiIntPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatRoffdiag = integrateRefEdgeTetraPhiIntPhiExtPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatSdiag = integrateRefEdgeTetraPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatSoffdiag = integrateRefEdgeTetraPhiIntPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad);

%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH = assembleMatElemDphiPhi(problemData.g, hatH);
problemData.globQ = assembleMatEdgeTetraPhiPhiNu(problemData.g, problemData.g.markE0Tint, hatSdiag, hatSoffdiag);
problemData.globQN = assembleMatEdgeTetraPhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, hatSdiag);
problemData.globS = problemData.eta * assembleMatEdgeTetraPhiPhi(problemData.g, problemData.g.markE0Tint, hatSdiag, hatSoffdiag, ones(problemData.g.numT, 4));
problemData.globSD = problemData.eta * assembleMatEdgeTetraPhiIntPhiInt(problemData.g, problemData.g.markE0TbdrD | problemData.g.markE0TbdrCoupling, hatSdiag, ones(problemData.g.numT, 4));

if ~problemData.isStationary
  problemData.sysW = [ sparse(2 * problemData.g.numT * problemData.N, 3 * problemData.g.numT * problemData.N) ; ...
                       sparse(problemData.g.numT * problemData.N, 2 * problemData.g.numT * problemData.N), problemData.globM ];
end % if

%% Empty vectors for coupled problem
problemData.globJcouple = { sparse(problemData.g.numT * problemData.N, 1), sparse(problemData.g.numT * problemData.N, 1) };
problemData.globKcouple = sparse(problemData.g.numT * problemData.N, 1);
end % function