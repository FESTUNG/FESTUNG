function problemData = outputStep(problemData, nStep)
%% visualization
if problemData.isVisSol && mod(nStep, problemData.outputFrequency) == 0
  cLagrange = projectDataDisc2DataLagr(problemData.cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, ...
                    nStep, problemData.outputTypes);
end % if
end % function

