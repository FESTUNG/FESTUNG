function problemData = outputStep(problemData, nStep)
%% visualization
varName = {};
dataLagr = {};
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species} && mod(nStep, problemData.outputFrequency{species}) == 0
    varName = [varName, {['u_' num2str(species)]}]; %#ok<AGROW>
    dataLagr = [dataLagr, {projectDataDisc2DataLagr(problemData.cDisc{species})}]; %#ok<AGROW>
  end % if
end % for
if ~isempty(dataLagr)
  visualizeDataLagr(problemData.g, dataLagr, varName, problemData.outputBasename, nStep, problemData.outputTypes)
end % if
end % function

