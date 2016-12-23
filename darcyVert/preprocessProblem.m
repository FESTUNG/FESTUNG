function problemData = preprocessProblem(problemData)

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
if problemData.isVisGrid, visualizeGridTrap(problemData.g); end

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Additional mesh data
problemData.g.markE0Tint = problemData.g.idE0T == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrN = problemData.g.idE0T == 2 | problemData.g.idE0T == 4; % [K x 3] mark local edges on the Neumann boundary
problemData.g.markE0TbdrD = ~(problemData.g.markE0Tint | problemData.g.markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d trapezoidals.\n', problemData.p, problemData.N, problemData.g.numT)

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuadTrap(problemData.p, problemData.qOrd);

%% Computation of matrices on the reference element.
problemData.hatM = integrateRefElemTrapPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatG = integrateRefElemTrapDphiPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatH = integrateRefElemTrapDphiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatRdiag = integrateRefEdgeTrapPhiIntPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatRoffdiag = integrateRefEdgeTrapPhiIntPhiExtPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatSdiag = integrateRefEdgeTrapPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
hatSoffdiag = integrateRefEdgeTrapPhiIntPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad);

%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH = assembleMatElemDphiPhi(problemData.g, hatH);
problemData.globQ = assembleMatEdgeTrapPhiPhiNu(problemData.g, problemData.g.markE0Tint, hatSdiag, hatSoffdiag);
problemData.globQN = assembleMatEdgeTrapPhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, hatSdiag);
problemData.globS = problemData.eta * assembleMatEdgeTrapPhiPhi(problemData.g, problemData.g.markE0Tint, hatSdiag, hatSoffdiag);
problemData.globSD = problemData.eta * assembleMatEdgeTrapPhiIntPhiInt(problemData.g, problemData.g.markE0TbdrD, hatSdiag);

problemData.sysW = [ sparse(2 * problemData.g.numT * problemData.N, 3 * problemData.g.numT * problemData.N) ; ...
                     sparse(problemData.g.numT * problemData.N, 2 * problemData.g.numT * problemData.N), problemData.globM ];
end % function