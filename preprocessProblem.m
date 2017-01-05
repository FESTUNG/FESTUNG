function problemData = preprocessProblem(problemData)

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
if problemData.isVisGrid, execin('darcyVert/visualizeGridTrap', problemData.g); end

%% Additional mesh data
problemData.g.markE0Tint = problemData.g.idE0T == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0Tbdr = ~problemData.g.markE0Tint; % all boundaries
problemData.g.markE0TbdrL = problemData.g.idE0T == 2 | problemData.g.idE0T == 4; % land boundaries
problemData.g.markE0TbdrF = problemData.g.idE0T == 3; % free boundary
problemData.g.markE0TbdrB = problemData.g.idE0T == 1; % bottom boundary

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs on trapezoidals
problemData.barN = problemData.p + 1;  % number of local DOFs on intervals
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d trapezoidals.\n', problemData.p, problemData.N, problemData.g.numT)

%% Lookup table for basis function.
problemData.basesOnQuad1D = computeBasesOnQuad1D(problemData.p, problemData.qOrd);
problemData.basesOnQuad2D = execin('darcyVert/computeBasesOnQuadTrap', problemData.p, problemData.qOrd);

%% Computation of matrices on the reference element.
problemData.hatM = execin('darcyVert/integrateRefElemTrapPhiPhi', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatG = execin('darcyVert/integrateRefElemTrapDphiPhiPhi', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
hatH = execin('darcyVert/integrateRefElemTrapDphiPhi', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
hatQdiag = execin('darcyVert/integrateRefEdgeTrapPhiIntPhiInt', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
hatQoffdiag = execin('darcyVert/integrateRefEdgeTrapPhiIntPhiExt', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatRdiag = execin('darcyVert/integrateRefEdgeTrapPhiIntPhiIntPhiInt', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatRoffdiag = execin('darcyVert/integrateRefEdgeTrapPhiIntPhiExtPhiExt', problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
hatSdiag = integrateRefEdgeTrapPhiIntPerQuad(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

problemData.barHatM = integrateRefElem1DPhiPhi(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barHatG = integrateRefElem1DDphiPhiPhi(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

barHatH = integrateRefElemTrapDphiPhi1D([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
barHatQdiag = integrateRefEdgeTrapPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
barHatQoffdiag = integrateRefEdgeTrapPhiIntPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.tildeHatPdiag = integrateRefEdgeTrapPhiIntPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.tildeHatPoffdiag = integrateRefEdgeTrapPhiIntPhiExtPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);

%% Assembly of time-independent global matrices.
% TODO: Adapt free surface
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH = assembleMatElemDphiPhi(problemData.g, hatH);
problemData.globQ = execin('darcyVert/assembleMatEdgeTrapPhiPhiNu', problemData.g, problemData.g.markE0Tint, hatQdiag, hatQoffdiag);

globQAvg = execin('darcyVert/assembleMatEdgeTrapPhiPhiNu', problemData.g, problemData.g.markE0Tint, hatQdiag, hatQoffdiag, 1:2);
problemData.globQavg = globQAvg{1};
problemData.globQup = assembleMatEdgeTrapPhiPhiNuBottomUp(problemData.g, problemData.g.markE0Tint | problemData.g.markE0TbdrF, hatQdiag, hatQoffdiag);

problemData.globS = assembleMatEdgeTrapPhiPerQuad(problemData.g, hatSdiag);

problemData.barGlobH = assembleMatElemDphiPhi1D(problemData.g, barHatH);
problemData.barGlobH = cellfun(@(c) problemData.gConst * c, problemData.barGlobH, 'UniformOutput', false);
problemData.barGlobQ = assembleMatEdgeTrapPhiPhi1DNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, barHatQdiag, barHatQoffdiag);
problemData.barGlobQbdr = assembleMatEdgeTrapPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tbdr, barHatQdiag);

problemData.barGlobM = assembleMatElemPhiPhi(problemData.g.g1D, problemData.barHatM);
end % function