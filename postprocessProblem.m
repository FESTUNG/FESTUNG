function problemData = postprocessProblem(problemData)
%% Visualization and error evaluation
varName = {};
dataLagr = {};
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species}
    varName = [varName, {['u_' num2str(species)]}]; %#ok<AGROW>
    dataLagr = [dataLagr, {projectDataDisc2DataLagr(problemData.cDisc{species})}]; %#ok<AGROW>
  end % if
  
  fprintf('Species %d L2 error w.r.t. the initial condition: %g\n', species, ...
    computeL2Error(problemData.g, problemData.cDisc{species}, problemData.c0Cont{species}, ...
                   2*problemData.p, problemData.basesOnQuad) );
end % for
if ~isempty(dataLagr)
  visualizeDataLagr(problemData.g, dataLagr, varName, problemData.outputBasename, problemData.numSteps, problemData.outputTypes)
end % if
end % function

