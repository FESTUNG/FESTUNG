function problemData = initializeProblem(problemData)
problemData.isFinished = false;

% Vector of unknowns (H, U, W) for each Runge-Kutta stage
problemData.cDiscRK = cell(length(rungeKuttaExplicit(problemData.ordRK, 0, 0)), 3);

%% Initial height.
problemData.cDiscRK{end, 1} = projectFuncCont2DataDisc1D(problemData.g.g1D, problemData.h0Cont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D);

%% Mesh adaptivity and assembly of time-independent global matrices.
problemData = problemData.fn_adaptFreeSurface(problemData, true);

%% Computation of bathymetry gradient.
dZbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 2);
dXbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 1) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 1);
dXzBot = problemData.g.g1D.markT2DT * ( problemData.gConst * (dZbot1D ./ dXbot1D) );
problemData.globLzBot = kron(dXzBot, eye(problemData.N, 1));

%% Initial velocities.
problemData.cDiscRK{end, 2} = projectFuncCont2DataDiscTetra(problemData.g, problemData.u10Cont, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
problemData.cDiscRK{end, 3} = zeros(size(problemData.cDiscRK{end, 2}));

%% Initial error computation.
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

fprintf('L2 errors of h, u1 w.r.t. the initial condition: %g, %g\n', ...
  computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{end, 1}, problemData.h0Cont, problemData.qOrd+1, problemData.basesOnQuad1D), ...
  computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 2}, problemData.u10Cont, problemData.qOrd+1, problemData.basesOnQuad2D));
end % function