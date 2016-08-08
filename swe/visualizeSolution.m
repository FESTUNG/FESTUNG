function visualizeSolution(pd, nStep)
if mod(nStep, pd.outputFrequency) == 0 || pd.isFinished
  nOutput = ceil(nStep / pd.outputFrequency);
  varName = {};
  dataLagr = {};
  
  %% Depth and elevation
  hDisc = pd.cDisc(:,:,1) - pd.zbDisc;

  % Evaluate water depth (h)
  if any(ismember(pd.outputList, 'h'))
    varName = [ varName, {'h'} ];
    dataLagr = [ dataLagr, {projectDataDisc2DataLagr(hDisc)} ];
  end % if

  % Evaluate water elevation (xi)
  if any(ismember(pd.outputList, 'xi')) || isfield(pd, 'stationElev')
    varName = [ varName, {'xi'} ];
    dataLagr = [ dataLagr, {projectDataDisc2DataLagr(pd.cDisc(:,:,1))} ];
  end % if

  % Save elevation station values
  if isfield(pd, 'stationElev')
    for n = 1 : length(pd.stationElev)
      dataStationV0T = dataLagr{end}(pd.stationElev{n}(:,1),:); % Extract values in vertices of relevant triangles
      pd.dataElev{n} = [ pd.dataElev{n} ; ...     % Append mean of barycentric weighted values
                         mean(sum(pd.stationElev{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if

  %% Momentum 
  if any(ismember(pd.outputList, {'uH', 'vH', 'momentum'}))
    varName = [ varName, {'uH', 'vH'} ];
    dataLagr = [ dataLagr, {projectDataDisc2DataLagr(pd.cDisc(:,:,2)), projectDataDisc2DataLagr(pd.cDisc(:,:,3))} ];
  end % if
  
  %% Velocity
  % Evaluate velocity
  if any(ismember(pd.outputList, {'u', 'v', 'velocity'})) || isfield(pd, 'stationVel')
    varName = [ varName, {'u', 'v'} ];
    dataQ0T = (pd.cDisc(:,:,2) * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.') ./ (hDisc * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.');
    dataDisc = pd.swe_projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
    dataLagr = [ dataLagr, {projectDataDisc2DataLagr(dataDisc)} ];
    
    dataQ0T = (pd.cDisc(:,:,3) * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.') ./ (hDisc * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.');
    dataDisc = pd.swe_projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
    dataLagr = [ dataLagr, {projectDataDisc2DataLagr(dataDisc)} ];
  end % if

  % Save velocity station values
  if isfield(pd, 'stationVel')
    for n = 1 : length(pd.stationVel)
      dataStationV0T = dataLagr{end-1}(pd.stationVel{n}(:,1),:); % Extract values in vertices of relevant triangles
      pd.dataVel{n,1} = [ pd.dataVel{n,1} ; ...     % Append mean of barycentric weighted values
                          mean(sum(pd.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
                        
      dataStationV0T = dataLagr{end}(pd.stationVel{n}(:,1),:); % Extract values in vertices of relevant triangles
      pd.dataVel{n,2} = [ pd.dataVel{n,1} ; ...     % Append mean of barycentric weighted values
                          mean(sum(pd.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
    end % for
  end % if
  
  %% Write visualization output
  outputMask = ismember(varName, pd.outputList);
  if any(outputMask)
    visualizeDataLagr(pd.g, dataLagr(outputMask), varName(outputMask), ['output/' pd.name], nOutput, pd.outputTypes);
  end % if
end % if
end % function

