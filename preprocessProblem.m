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

problemData.g.g1D.markV0Tint = problemData.g.g1D.idV0T == 0;
problemData.g.g1D.markV0Tbdr = ~problemData.g.g1D.markV0Tint;

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
problemData.barHatG = integrateRefElem1DDphiPhiPhiPerQuad(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
barHatSdiag = integrateRefEdgeTrapPhi1DIntPerQuad(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barHatPdiag = integrateRefEdge1DPhiIntPhiIntPhiInt(problemData.barN, problemData.basesOnQuad1D);
problemData.barHatPoffdiag = integrateRefEdge1DPhiIntPhiExtPhiExt(problemData.barN, problemData.basesOnQuad1D);

tildeHatH = integrateRefElemTrapDphiPhi1D([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
tildeHatQdiag = integrateRefEdgeTrapPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
tildeHatQoffdiag = integrateRefEdgeTrapPhiIntPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
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
problemData.barGlobS = assembleMatEdgeTrapPhi1DPerQuad(problemData.g, barHatSdiag);

problemData.tildeGlobH = assembleMatElemDphiPhi1D(problemData.g, tildeHatH);
problemData.tildeGlobQ = assembleMatEdgeTrapPhiPhi1DNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, tildeHatQdiag, tildeHatQoffdiag);
problemData.tildeGlobQbdr = assembleMatEdgeTrapPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tbdr, tildeHatQdiag);
for m = 1 : 2
  problemData.tildeGlobH = cellfun(@(c) problemData.gConst * c, problemData.tildeGlobH, 'UniformOutput', false);
  problemData.tildeGlobQ = cellfun(@(c) problemData.gConst * c, problemData.tildeGlobQ, 'UniformOutput', false);
  problemData.tildeGlobQbdr = cellfun(@(c) problemData.gConst * c, problemData.tildeGlobQbdr, 'UniformOutput', false);
end % for m

problemData.barGlobM = assembleMatElemPhiPhi(problemData.g.g1D, problemData.barHatM);
end % function