function visualizeSolution(pd, nStep)
if mod(nStep, pd.outputFrequency) == 0 || pd.isFinished
  nOutput = ceil(nStep / pd.outputFrequency);

  %% Depth and elevation
  hDisc = pd.cDisc(:,:,1) - pd.zbDisc;

  % Visualize water depth (h)
  if any(ismember(pd.outputList, 'h'))
    dataLagr = projectDataDisc2DataLagr(hDisc);
    visualizeDataLagr(pd.g, dataLagr, 'h_h', ['output/' pd.name '_h'], nOutput, pd.outputTypes);
  end % if

  % Evaluate water elevation (xi)
  if any(ismember(pd.outputList, 'xi')) || isfield(pd, 'stationElev')
    dataLagr = projectDataDisc2DataLagr(pd.cDisc(:,:,1));
  end % if

  % Visualize water elevation (xi)
  if any(ismember(pd.outputList, 'xi'))
    visualizeDataLagr(pd.g, dataLagr, 'xi_h', ['output/' pd.name '_xi'], nOutput, pd.outputTypes);
  end % if

  % Save elevation station values
  if isfield(pd, 'stationElev')
    for n = 1 : length(pd.stationElev)
      dataStationV0T = dataLagr(pd.stationElev{n}(:,1),:); % Extract values in vertices of relevant triangles
      pd.dataElev{n} = [ pd.dataElev{n} ; ...     % Append mean of barycentric weighted values
                         mean(sum(pd.stationElev{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if

  %% Momentum and velocity (first component)

  % Visualize primary variable uH
  if any(ismember(pd.outputList, 'uH'))
    dataLagr = projectDataDisc2DataLagr(pd.cDisc(:,:,2));
    visualizeDataLagr(pd.g, dataLagr, 'uH_h', ['output/' pd.name '_uH'], nOutput, pd.outputTypes);
  end % if

  % Evaluate x-velocity (u)
  if any(ismember(pd.outputList, 'u')) || isfield(pd, 'stationVel')
    dataQ0T = (pd.cDisc(:,:,2) * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.') ./ (hDisc * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.');
    dataDisc = pd.swe_projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
  end % if

  % Visualize x-velocity (u)
  if any(ismember(pd.outputList, 'u'))
    visualizeDataLagr(pd.g, dataLagr, 'u_h', ['output/' pd.name '_u'], nOutput, pd.outputTypes);
  end % if

  % Save x-velocity station values
  if isfield(pd, 'stationVel')
    for n = 1 : length(pd.stationVel)
      dataStationV0T = dataLagr(pd.stationVel{n}(:,1),:); % Extract values in vertices of relevant triangles
      pd.dataVel{n,1} = [ pd.dataVel{n,1} ; ...     % Append mean of barycentric weighted values
                          mean(sum(pd.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if

  %% Momentum and velocity (second component)

  % Visualize primary variable vH
  if any(ismember(pd.outputList, 'vH'))
    dataLagr = projectDataDisc2DataLagr(pd.cDisc(:,:,3));
    visualizeDataLagr(pd.g, dataLagr, 'vH_h', ['output/' pd.name '_vH'], nOutput, pd.outputTypes);
  end % if

  % Evaluate y-velocity (v)
  if any(ismember(pd.outputList, 'v')) || isfield(pd, 'stationVel')
    dataQ0T = (pd.cDisc(:,:,3) * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.') ./ (hDisc * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.');
    dataDisc = pd.swe_projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
  end % if

  % Visualize y-velocity (v)
  if any(ismember(pd.outputList, 'v'))
    visualizeDataLagr(pd.g, dataLagr, 'v_h', ['output/' pd.name '_v'], nOutput, pd.outputTypes);
  end % if

  % Save y-velocity station values
  if isfield(pd, 'stationVel')
    for n = 1 : length(pd.stationVel)
      dataStationV0T = dataLagr(pd.stationVel{n}(:,1),:); % Extract values in vertices of relevant triangles
      pd.dataVel{n,2} = [ pd.dataVel{n,2} ; ...     % Append mean of barycentric weighted values
                          mean(sum(pd.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if
end % if
end % function

