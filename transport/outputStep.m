function problemData = outputStep(problemData, nStep)
%% visualization
hQ0T = (problemData.hDisc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.');
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species} && mod(nStep, problemData.outputFrequency{species}) == 0
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['cH_' num2str(species) '_h'], ...
                      [problemData.outputBasename{species} '_H'], ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.') ./ hQ0T;
    dataDisc = problemData.swe_projectDataQ0T2DataDisc(dataQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
    cLagrange = projectDataDisc2DataLagr(dataDisc);
    visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                      problemData.outputBasename{species}, ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
  end % if
end % for
if ~isempty(dataLagr)
  visualizeDataLagr(problemData.g, dataLagr, varName, problemData.outputBasename, nStep, problemData.outputTypes)
end % if
end % function
