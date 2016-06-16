function problemData = postprocessProblem(problemData)
for species = 1:problemData.numSpecies
  %% Visualization
  if problemData.isVisSol{species}
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['u_' num2str(species) '_h'], problemData.outputBasename{species}, ...
                      problemData.numSteps, problemData.outputTypes{species});
  end % if
  %% Error evaluation
  fprintf('L2 error w.r.t. the initial condition: %g\n', ...
    computeL2Error(problemData.g, problemData.cDisc{species}, problemData.c0Cont{species}, 2*problemData.p, problemData.basesOnQuad));
  fprintf('norm(cDisc, 1) = %g\n', norm(problemData.cDisc{species}(:), 1));
  fprintf('norm(cDisc, 2) = %g\n', norm(problemData.cDisc{species}(:), 2));
  fprintf('norm(cDisc, inf) = %g\n', norm(problemData.cDisc{species}(:), inf));
end % for
end % function

