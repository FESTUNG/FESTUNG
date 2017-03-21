function problemData = initializeProblem(problemData)
problemData.isFinished = false;
problemData.cDisc = cell(1,3);  % Vector of unknowns (H, U, W)

%% Initial height.
h0Cont = @(x1) problemData.hCont(0,x1);
problemData.cDisc{1} = projectFuncCont2DataDisc1D(problemData.g.g1D, h0Cont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D);

%% Mesh adaptivity and assembly of time-independent global matrices.
problemData = adaptMesh(problemData, true);

%% Initial velocity.
u10Cont = @(x1,x2) problemData.u1Cont(0,x1,x2);
u20Cont = @(x1,x2) problemData.u2Cont(0,x1,x2);
problemData.cDisc{2} = problemData.fn_projectFuncCont2DataDiscTrap(problemData.g, u10Cont, problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
problemData.cDisc{3} = problemData.fn_projectFuncCont2DataDiscTrap(problemData.g, u20Cont, problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);

%% Error computation and visualization of inital condition.
if problemData.isVisGrid, execin('darcyVert/visualizeGridTrap', problemData.g); end

fprintf('L2 errors of cDisc w.r.t. the initial condition: %g, %g, %g\n', ...
  computeL2Error1D(problemData.g.g1D, problemData.cDisc{1}, h0Cont, problemData.qOrd+1, problemData.basesOnQuad1D), ...
  execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDisc{2}, u10Cont, problemData.qOrd+1, problemData.basesOnQuad2D), ...
  execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDisc{3}, u20Cont, problemData.qOrd+1, problemData.basesOnQuad2D));

if problemData.isVisSol
  cLagr = cellfun(@(c) execin('darcyVert/projectDataDisc2DataLagrTrap', c), problemData.cDisc(2:3), 'UniformOutput', false);
  problemData.fn_visualizeDataLagrTrap(problemData.g, cLagr, {'u1', 'u2'}, problemData.outputBasename, 0, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
end % if
end % function