function problemData = initializeProblem(problemData)
problemData.isFinished = false;

% Vector of unknowns (H, U, W) for each Runge-Kutta stage
problemData.cDiscRK = cell(length(rungeKuttaSSP(problemData.ordRK, 0, 0)), 3);

%% Initial height.
problemData.cDiscRK{1, 1} = projectFuncCont2DataDisc1D(problemData.g.g1D, problemData.h0Cont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D);

%% Mesh adaptivity and assembly of time-independent global matrices.
problemData = adaptMesh(problemData, true);

%% Computation of bathymetry gradient.
dZbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 2);
dXbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 1) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 1);
dXzBot = problemData.g.g1D.markT2DT * ( problemData.gConst * (dZbot1D ./ dXbot1D) );
problemData.globLzBot = kron(dXzBot, eye(problemData.N, 1));

%% Initial velocities.
problemData.cDiscRK{1, 2} = projectFuncCont2DataDiscTetra(problemData.g, problemData.u10Cont, problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
problemData.cDiscRK{1, 3} = zeros(size(problemData.cDiscRK{1, 2}));

%% Error computation and visualization of inital condition.
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

fprintf('L2 errors of h, u1 w.r.t. the initial condition: %g, %g\n', ...
  computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{1, 1}, problemData.h0Cont, problemData.qOrd+1, problemData.basesOnQuad1D), ...
  computeL2ErrorTetra(problemData.g, problemData.cDiscRK{1, 2}, problemData.u10Cont, problemData.qOrd+1, problemData.basesOnQuad2D));

if problemData.isVisSol
  cLagr = cellfun(@(c) projectDataDisc2DataLagrTensorProduct(c), problemData.cDiscRK(1, 2:3), 'UniformOutput', false);
  visualizeDataLagrTetra(problemData.g, cLagr, {'u1', 'u2'}, problemData.outputBasename, 0, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
end % if
end % function