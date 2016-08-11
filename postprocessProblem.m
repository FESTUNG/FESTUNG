function problemData = postprocessProblem(problemData)
problemData.errors = zeros(problemData.numSpecies,1);
for species = 1:problemData.numSpecies
%   %% Visualization
%   if problemData.isVisSol{species}
%     cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
%     visualizeDataLagr(problemData.g, cLagrange, ['cH_' num2str(species) '_h'], problemData.outputBasename{species}, ...
%                       problemData.numSteps, problemData.outputTypes{species});
%   end % if
  %% Error evaluation
%   fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  problemData.errors(species) = computeL2Error(problemData.g, problemData.cDisc{species}, @(x1,x2) problemData.solCont{species}(problemData.tEnd,x1,x2), 2*problemData.p, problemData.basesOnQuad);%);
%   fprintf('norm(cDisc, 1) = %g\n', norm(problemData.cDisc{species}(:), 1));
%   fprintf('norm(cDisc, 2) = %g\n', norm(problemData.cDisc{species}(:), 2));
%   fprintf('norm(cDisc, inf) = %g\n', norm(problemData.cDisc{species}(:), inf));
end % for
if ~isempty(dataLagr)
  visualizeDataLagr(problemData.g, dataLagr, varName, problemData.outputBasename, problemData.numSteps, problemData.outputTypes)
end % if
end % function

