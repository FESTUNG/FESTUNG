function problemData = outputStep(problemData, nStep)
%% visualization
varName = {};
dataLagr = {};
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species} && mod(nStep, problemData.outputFrequency{species}) == 0
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['cH_' num2str(species) '_h'], ...
                      problemData.outputBasename{species}, ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
  end % if
end % for
if ~isempty(dataLagr)
  visualizeDataLagr(problemData.g, dataLagr, varName, problemData.outputBasename, nStep, problemData.outputTypes)
end % if
end % function
