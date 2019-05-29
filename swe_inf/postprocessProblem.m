
function pd = postprocessProblem(pd)
%% Close waitbar.
if pd.isWaitbar
  close(pd.waitbar)
end % if


%% Compute error if analytical solution available.
p = pd.p;
if pd.isSolutionAvail
  % Continuous solution
  t = pd.t0 + pd.numSteps * pd.dt;
  
  % Error in water height (H)
  [errHInf, erruHInf, errvHInf] = computeErrorProjection(pd.gInf, @(x,y) pd.HOSAlg(x,y,t), @(x,y) pd.uCont(x,y,t), @(x,y) pd.vCont(x,y,t), pd.cDiscInf, pd.p, pd.pInf, pd.beta);
  fprintf('L2-Error H: %g\n', errHInf);
  fprintf('L2-Error uH: %g\n', erruHInf);
  fprintf('L2-Error vH: %g\n', errvHInf);
  
end % if
end

