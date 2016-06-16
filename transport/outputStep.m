function problemData = outputStep(problemData, nStep)
%% visualization
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species} && mod(nStep, problemData.outputFrequency{species}) == 0
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['u_' num2str(species) '_h'], problemData.outputBasename{species}, ...
                      nStep, problemData.outputTypes{species});
  end % if
end % for
end % function

