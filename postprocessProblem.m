function problemData = postprocessProblem(problemData)
problemData.errors = zeros(problemData.numSpecies,2);

hQ0T = (problemData.hDisc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.');

if problemData.isSolutionAvailable
  for species = 1:problemData.numSpecies
    %% Error evaluation
    problemData.errors(species,1) = computeL2Error(problemData.g, problemData.cDisc{species}, @(x1,x2) problemData.solCont{species}(problemData.tEnd,x1,x2) .* problemData.hCont(x1,x2,problemData.tEnd), 2*problemData.p, problemData.basesOnQuad);
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.') ./ hQ0T;
    dataDisc = problemData.swe_projectDataQ0T2DataDisc(dataQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
    problemData.errors(species,2) = computeL2Error(problemData.g, dataDisc, @(x1,x2) problemData.solCont{species}(problemData.tEnd,x1,x2), 2*problemData.p, problemData.basesOnQuad);
  end % for
end % if
end % function

