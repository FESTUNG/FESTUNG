function visualizeSolution(problemData, nStep)
p = problemData.p;

if p <= 2 && mod(nStep, problemData.outputFrequency) == 0
  nOutput = nStep / problemData.outputFrequency;

  %% Depth and elevation

  % Visualize water depth (h)
  if any(ismember(problemData.outputList, 'h'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,1));
    visualizeDataLagr(problemData.g, dataLagr, 'h_h', ['output/' problemData.name '_h'], nOutput, problemData.outputTypes);
  end % if

  % Evaluate water elevation (xi)
  if any(ismember(problemData.outputList, 'xi')) || ~isempty(problemData.stationElev)
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,1) + problemData.zbDisc);
  end % if

  % Visualize water elevation (xi)
  if any(ismember(problemData.outputList, 'xi'))
    visualizeDataLagr(problemData.g, dataLagr, 'xi_h', ['output/' problemData.name '_xi'], nOutput, problemData.outputTypes);
  end % if

  % Save elevation station values
  if ~isempty(problemData.stationElev)
    for n = 1 : length(problemData.stationElev)
      dataStationV0T = dataLagr(problemData.stationElev{n}(:,1),:); % Extract values in vertices of relevant triangles
      problemData.dataElev{n} = [ problemData.dataElev{n} ; ...     % Append mean of barycentric weighted values
                                  mean(sum(problemData.stationElev{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if

  %% Momentum and velocity (first component)

  % Visualize primary variable uH
  if any(ismember(problemData.outputList, 'uH'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,2));
    visualizeDataLagr(problemData.g, dataLagr, 'uH_h', ['output/' problemData.name '_uH'], nOutput, problemData.outputTypes);
  end % if

  % Evaluate x-velocity (u)
  if any(ismember(problemData.outputList, 'u')) || ~isempty(problemData.stationVel)
    dataDisc = projectQuotientDisc2Disc(problemData.cDisc(:,:,2), problemData.cDisc(:,:,1), 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
  end % if

  % Visualize x-velocity (u)
  if any(ismember(problemData.outputList, 'u'))
    visualizeDataLagr(problemData.g, dataLagr, 'u_h', ['output/' problemData.name '_u'], nOutput, problemData.outputTypes);
  end % if

  % Save x-velocity station values
  if ~isempty(problemData.stationVel)
    for n = 1 : length(problemData.stationVel)
      dataStationV0T = dataLagr(problemData.stationVel{n}(:,1),:); % Extract values in vertices of relevant triangles
      problemData.dataVel{n,1} = [ problemData.dataVel{n,1} ; ...     % Append mean of barycentric weighted values
                                   mean(sum(problemData.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if

  %% Momentum and velocity (second component)

  % Visualize primary variable vH
  if any(ismember(problemData.outputList, 'vH'))
    dataLagr = projectDataDisc2DataLagr(problemData.cDisc(:,:,3));
    visualizeDataLagr(problemData.g, dataLagr, 'vH_h', ['output/' problemData.name '_vH'], nOutput, problemData.outputTypes);
  end % if

  % Evaluate y-velocity (v)
  if any(ismember(problemData.outputList, 'v')) || ~isempty(problemData.stationVel)
    dataDisc = projectQuotientDisc2Disc(problemData.cDisc(:,:,3), problemData.cDisc(:,:,1), 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
  end % if

  % Visualize y-velocity (v)
  if any(ismember(problemData.outputList, 'v'))
    visualizeDataLagr(problemData.g, dataLagr, 'v_h', ['output/' problemData.name '_v'], nOutput, problemData.outputTypes);
  end % if

  % Save y-velocity station values
  if ~isempty(problemData.stationVel)
    for n = 1 : length(problemData.stationVel)
      dataStationV0T = dataLagr(problemData.stationVel{n}(:,1),:); % Extract values in vertices of relevant triangles
      problemData.dataVel{n,2} = [ problemData.dataVel{n,2} ; ...     % Append mean of barycentric weighted values
                                   mean(sum(problemData.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if
end

