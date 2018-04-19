function problemData = initializeProblem(problemData)
problemData.isFinished = false;
K = problemData.g.numT;
N = problemData.N;

%% Initial data
if problemData.isHotstart
  hotstart = load(problemData.hotstartFile);
  assert(all(isfield(hotstart, {'t', 'hDisc', 'q1Disc', 'q2Disc'})), 'Hotstart file must contain t, hDisc, q1Disc, q2Disc')
  assert(all(problemData.g.numT == [size(hotstart.hDisc, 1), size(hotstart.q1Disc, 1), size(hotstart.q2Disc, 1)]), ...
    'Number of elements in hotstart file does not match!')
  fprintf('Loaded hotstart data from "%s" at time level t=%g.\n', problemData.hotstartFile, hotstart.t);
  
  hDisc = zeros(K, N);
  q1Disc = zeros(K, N);
  q2Disc = zeros(K, N);
  
  hotstartN = min(N, size(hotstart.hDisc, 2));
  hDisc(:, 1:hotstartN) = hotstart.hDisc(:, 1:hotstartN);
  q1Disc(:, 1:hotstartN) = hotstart.q1Disc(:, 1:hotstartN);
  q2Disc(:, 1:hotstartN) = hotstart.q2Disc(:, 1:hotstartN);
else
  if isfield(problemData, 'hCont')
    h0Cont = @(x1,x2) problemData.hCont(0,x1,x2);
  elseif isfield(problemData, 'h0Cont')
    h0Cont = problemData.h0Cont;
  else
    warning('No initial data for h given. Initializing with zeros.');
    h0Cont = @(x1,x2) zeros(size(x1));
  end % if

  if isfield(problemData, 'q1Cont')
    q10Cont = @(x1,x2) problemData.q1Cont(0,x1,x2);
  elseif isfield(problemData, 'q10Cont')
    q10Cont = problemData.q10Cont;
  else
    warning('No initial data for q1 given. Initializing with zeros.');
    q10Cont = @(x1,x2) zeros(size(x1));
  end % if

  if isfield(problemData, 'q2Cont')
    q20Cont = @(x1,x2) problemData.q2Cont(0,x1,x2);
  elseif isfield(problemData, 'q20Cont')
    q20Cont = problemData.q20Cont;
  else
    warning('No initial data for q2 given. Initializing with zeros.');
    q20Cont = @(x1,x2) zeros(size(x1));
  end % if

  hDisc = projectFuncCont2DataDiscTetra(problemData.g, h0Cont, problemData.qOrd, ...
                                        problemData.globM, problemData.basesOnQuad);
  q1Disc = projectFuncCont2DataDiscTetra(problemData.g, q10Cont, problemData.qOrd, ...
                                         problemData.globM, problemData.basesOnQuad);
  q2Disc = projectFuncCont2DataDiscTetra(problemData.g, q20Cont, problemData.qOrd, ...
                                         problemData.globM, problemData.basesOnQuad);
                                       
  fprintf('L2 error w.r.t. the initial condition: %g\n', ...
    computeL2ErrorTetra(problemData.g, hDisc, h0Cont, problemData.qOrd+1, problemData.basesOnQuad));
end % if

problemData.sysY = [ reshape(q1Disc', K * N, 1) ; ...
                     reshape(q2Disc', K * N, 1) ; ...
                     reshape(hDisc', K * N, 1) ];

%% Visualization of inital condition.
if problemData.isVisSol
  cLagr = { projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)'), ...
            projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(1 : K*N), N, K)'), ...
            projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(K*N+1 : 2*K*N), N, K)') };
  visualizeDataLagrTetra(problemData.g, cLagr, {'h', 'q1', 'q2'}, problemData.outputBasename, 0, problemData.outputTypes, struct('q', {{'q1','q2'}}));
end % if

end % function