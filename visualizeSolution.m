function visualizeSolution(problemData, nStep)
p = problemData.p;

if p <= 2 && mod(nStep, problemData.outputFrequency) == 0
  nOutput = nStep / problemData.outputFrequency;
  
  % Visualize water depth (h)
  if any(ismember(problemData.outputList, 'h'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,1));
    visualizeDataLagr(problemData.g, dataLagr, 'h_h', ['output/' problemData.name '_h'], nOutput, problemData.outputTypes);
  end % if
  
  % Visualize water elevation (xi)
  if any(ismember(problemData.outputList, 'xi'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,1) + problemData.zbDisc);
    visualizeDataLagr(problemData.g, dataLagr, 'xi_h', ['output/' problemData.name '_xi'], nOutput, problemData.outputTypes);
  end % if
  
  % Visualize primary variable uH
  if any(ismember(problemData.outputList, 'uH'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,2));
    visualizeDataLagr(problemData.g, dataLagr, 'uH_h', ['output/' problemData.name '_uH'], nOutput, problemData.outputTypes);
  end % if
  
  % Visualize x-velocity (u)
  if any(ismember(problemData.outputList, 'u'))
    dataDisc = projectQuotientDisc2Disc(problemData.cDisc(:,:,2), problemData.cDisc(:,:,1), 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
    visualizeDataLagr(problemData.g, dataLagr, 'u_h', ['output/' problemData.name '_u'], nOutput, problemData.outputTypes);
  end % if
  
  % Visualize primary variable vH
  if any(ismember(problemData.outputList, 'vH'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,3));
    visualizeDataLagr(problemData.g, dataLagr, 'vH_h', ['output/' problemData.name '_vH'], nOutput, problemData.outputTypes);
  end % if
  
  % Visualize y-velocity (v)
  if any(ismember(problemData.outputList, 'v'))
    dataDisc = projectQuotientDisc2Disc(problemData.cDisc(:,:,3), problemData.cDisc(:,:,1), 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
    visualizeDataLagr(problemData.g, dataLagr, 'v_h', ['output/' problemData.name '_v'], nOutput, problemData.outputTypes);
  end % if
end

