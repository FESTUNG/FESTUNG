function problemData = postprocessProblem(problemData)
%% Visualization
if problemData.isVisSol
  cLagrange = projectDataDisc2DataLagr(problemData.cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, ...
                    problemData.numSteps, problemData.outputTypes);
end % if
%% Error evaluation
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
        computeL2Error(problemData.g, problemData.cDisc, problemData.c0Cont, 2*problemData.p));
fprintf('norm(cDisc, 1) = %g\n', norm(problemData.cDisc(:), 1));
fprintf('norm(cDisc, 2) = %g\n', norm(problemData.cDisc(:), 2));
fprintf('norm(cDisc, inf) = %g\n', norm(problemData.cDisc(:), inf));
end % function

