function problemData = initializeProblem(problemData)
problemData.isFinished = false;

% Vector of unknowns (H, U, W) for each Runge-Kutta stage
problemData.cDiscRK = cell(length(rungeKuttaExplicit(problemData.ordRK, 0, 0)), 3);

%% Load state
if problemData.isHotstart
  hotstart = load(problemData.hotstartFile);
  assert(all(isfield(hotstart, {'t', 'hDisc', 'u1Disc', 'u2Disc'})), 'Hotstart file must contain t, hDisc, u1Disc, u2Disc')
  assert(problemData.g.g1D.numT == size(hotstart.hDisc, 1), 'Number of 1D elements in hotstart file does not match!')
  assert(all(problemData.g.numT == [size(hotstart.u1Disc, 1), size(hotstart.u2Disc, 1)]), 'Number of 2D elements in hotstart file does not match!')
  fprintf('Loaded hotstart data from "%s" at time level t=%g.\n', problemData.hotstartFile, hotstart.t);
end % if

%% Initial height.
if problemData.isHotstart
  barN = min(problemData.barN, size(hotstart.hDisc, 2));
  problemData.cDiscRK{end, 1} = zeros(problemData.g.g1D.numT, problemData.barN);
  problemData.cDiscRK{end, 1}(:, 1:barN) = hotstart.hDisc(:, 1:barN);
else
  problemData.cDiscRK{end, 1} = projectFuncCont2DataDisc1D(problemData.g.g1D, problemData.h0Cont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D);
end % if

%% Mesh adaptivity and assembly of time-independent global matrices.
problemData = problemData.fn_adaptFreeSurface(problemData, true);

%% Computation of bathymetry gradient.
dZbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 2);
dXbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 1) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 1);
dXzBot = problemData.g.g1D.markT2DT * ( problemData.gConst * (dZbot1D ./ dXbot1D) );
problemData.globLzBot = kron(dXzBot, eye(problemData.N, 1));

%% Initial velocities.
if problemData.isHotstart
  N = min(problemData.N, size(hotstart.u1Disc, 2));
  problemData.cDiscRK{end, 2} = zeros(problemData.g.numT, problemData.N);
  problemData.cDiscRK{end, 2}(:, 1:N) = hotstart.u1Disc(:, 1:N);
  problemData.cDiscRK{end, 3} = zeros(problemData.g.numT, problemData.N);
  problemData.cDiscRK{end, 3}(:, 1:N) = hotstart.u2Disc(:, 1:N);
else
  problemData.cDiscRK{end, 2} = projectFuncCont2DataDiscTetra(problemData.g, problemData.u10Cont, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
  % problemData.cDiscRK{end, 3} = projectFuncCont2DataDiscTetra(problemData.g, @(x,z) problemData.u2Cont(0,x,z), problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
  problemData.cDiscRK{end, 3} = zeros(size(problemData.cDiscRK{end, 2}));
end % if

%% Exact solution
if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
  problemData.pdExact = problemData;
end % if

%% Initial error computation.
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

fprintf('L2 errors of h, u1 w.r.t. the initial condition: %g, %g\n', ...
  computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{end, 1}, problemData.h0Cont, problemData.qOrd+1, problemData.basesOnQuad1D), ...
  computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 2}, problemData.u10Cont, problemData.qOrd+1, problemData.basesOnQuad2D));
end % function