function problemData = preprocessProblem(problemData)

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
if problemData.isVisGrid, execin('darcyVert/visualizeGridTrap', problemData.g); end

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d trapezoidals.\n', problemData.p, problemData.N, problemData.g.numT)

%% Additional mesh data
problemData.g.markE0Tint = problemData.g.idE0T == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrL = problemData.g.idE0T == 2 | problemData.g.idE0T == 4; % land boundaries
problemData.g.markE0TbdrF = problemData.g.idE0T == 3; % free boundary
problemData.g.markE0TbdrB = problemData.g.idE0T == 1; % bottom boundary

%% Lookup table for basis function.
problemData.basesOnQuad = execin('darcyVert/computeBasesOnQuadTrap', problemData.p, problemData.qOrd);

%% Computation of matrices on the reference element.
problemData.hatM = execin('darcyVert/integrateRefElemTrapPhiPhi', problemData.N, problemData.qOrd, problemData.basesOnQuad);
% hatG = integrateRefElemTrapDphiPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatH = execin('darcyVert/integrateRefElemTrapDphiPhi', problemData.N, problemData.qOrd, problemData.basesOnQuad);
% problemData.hatRdiag = integrateRefEdgeTrapPhiIntPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
% problemData.hatRoffdiag = integrateRefEdgeTrapPhiIntPhiExtPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatSdiag = execin('darcyVert/integrateRefEdgeTrapPhiIntPhiInt', problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.hatSoffdiag = execin('darcyVert/integrateRefEdgeTrapPhiIntPhiExt', problemData.N, problemData.qOrd, problemData.basesOnQuad);

%% Assembly of time-independent global matrices.
% TODO: Adapt free surface
problemData.globM = execin('darcyVert/assembleMatElemTrapPhiPhi', problemData.g, problemData.hatM);
problemData.globH = execin('darcyVert/assembleMatElemTrapDphiPhi', problemData.g, problemData.hatH);
problemData.globQ = execin('darcyVert/assembleMatEdgeTrapPhiPhiNu', problemData.g, problemData.g.markE0Tint, problemData.hatSdiag, problemData.hatSoffdiag);

end % function