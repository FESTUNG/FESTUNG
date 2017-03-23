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

%% Initial horizontal velocity.
problemData.cDiscRK{1, 2} = problemData.fn_projectFuncCont2DataDiscTrap(problemData.g, problemData.u10Cont, problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);

%% Initial vertical velocity
% u1DCont = @(x1,x2) problemData.u1DCont(problemData.t0,x1,x2);
% u2DCont = @(x1,x2) problemData.u2DCont(problemData.t0,x1,x2);
% problemData.globJu = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
% problemData.globJw = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrW, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
% 
% [Q,~] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);
% heightV0T1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), [4 3], 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), [1 2], 2);
% 
% mapE0E = [2 1 4 3];
% u1Q0E0Tint = cell(4,1); % cDisc{2} in quad points of edges
% u1Q0E0TE0T = cell(4,1); % cDisc{2} of neighboring element in quad points of edges
% for n = 1 : 4
%   u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,n) * problemData.cDiscRK{1, 2}.', problemData.g.numT * numQuad1D, 1);
%   cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,mapE0E(n)) * problemData.cDiscRK{1, 2}.';
%   u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', problemData.g.numT * numQuad1D, 1);
% end % for nn
% 
% problemData.globKh = zeros(problemData.g.numT * problemData.N, 1);
% for n = 3 : 4
%   hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,5-n) + problemData.g.g1D.markV0TV0T{5-n} * problemData.hV0T1D(:,5-mapE0E(n)) );
%   hJmpE0T = problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,5-n) - problemData.g.g1D.markV0TV0T{5-n} * problemData.hV0T1D(:,5-mapE0E(n)) );
%   u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n});
%   lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
%   hJmpLambdaE0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D,1));
%     
%   problemData.globKh = problemData.globKh + problemData.globS{n} * hJmpLambdaE0T;
% end % for n
% 
% problemData.tildeGlobP = assembleMatEdgeTrapPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{1, 1}, heightV0T1D, problemData.g.markE0Tint, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);
% 
% problemData.cDiscRK{1, 3} = reshape( (problemData.globHQup) \ ( problemData.globJu{1} + problemData.globJw{2} + problemData.globKh + ...
%                                 (problemData.globHQavg + problemData.tildeGlobP) * reshape(problemData.cDiscRK{1, 2}.', [], 1) ), ...
%                              problemData.N, problemData.g.numT ).';

%% Error computation and visualization of inital condition.
if problemData.isVisGrid, execin('darcyVert/visualizeGridTrap', problemData.g); end

fprintf('L2 errors of h, u1 w.r.t. the initial condition: %g, %g\n', ...
  computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{1, 1}, problemData.h0Cont, problemData.qOrd+1, problemData.basesOnQuad1D), ...
  execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDiscRK{1, 2}, problemData.u10Cont, problemData.qOrd+1, problemData.basesOnQuad2D));

if problemData.isVisSol
  cLagr = cellfun(@(c) execin('darcyVert/projectDataDisc2DataLagrTrap', c), problemData.cDiscRK(1, 2:3), 'UniformOutput', false);
  problemData.fn_visualizeDataLagrTrap(problemData.g, cLagr, {'u1', 'u2'}, problemData.outputBasename, 0, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
end % if
end % function