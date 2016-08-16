function problemData = initializeProblem(problemData)
%% Initial data.
problemData.cDisc = cell(problemData.numSpecies,1);

hDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.hCont(0,x1,x2), 2*problemData.p+1,problemData.hatM, problemData.basesOnQuad);
hQ0T = (hDisc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.');

for species = 1:problemData.numSpecies
  problemData.cDisc{species} = projectFuncCont2DataDisc(problemData.g, problemData.cH0Cont{species}, 2*problemData.p+1, ...
                                                        problemData.hatM, problemData.basesOnQuad);
  if problemData.isSlopeLim{species}
    cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont{species}(0, x1, x2));
    problemData.cDisc{species} = applySlopeLimiterDisc(problemData.g, problemData.cDisc{species}, ...
                                                       problemData.g.markV0TbdrD, cDV0T, problemData.globM, ...
                                                       problemData.globMDiscTaylor, problemData.basesOnQuad, ...
                                                       problemData.typeSlopeLim{species});
  end % if
  
  % Initial error TODO maybe for non-integrated concentrations too?
  fprintf('L2 error w.r.t. the initial condition: %g\n', ...
    computeL2Error(problemData.g, problemData.cDisc{species}, problemData.cH0Cont{species}, 2*problemData.p, problemData.basesOnQuad));
  
  % Visualization of inital condition.
  if problemData.isVisSol{species}
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['cH_' num2str(species) '_h'], ... 
                      [problemData.outputBasename{species} '_H'], 0, problemData.outputTypes{species})
    
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.') ./ hQ0T;
    dataDisc = problemData.swe_projectDataQ0T2DataDisc(dataQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
    cLagrange = projectDataDisc2DataLagr(dataDisc);
    visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                      problemData.outputBasename{species}, 0, problemData.outputTypes{species});
  end % if
end % for
%% Initialize time stepping.
problemData.isFinished = false;
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', ...
  problemData.tEnd, problemData.tau, problemData.numSteps)
end % function