% Auxiliary function for visualization of 2D Shallow-Water Equations output.

%===============================================================================
%> @file
%>
%> @brief Auxiliary function for visualization of 2D Shallow-Water Equations 
%>        output.
%===============================================================================
%>
%> @brief Auxiliary function for visualization of 2D Shallow-Water Equations 
%>        output.
%>
%> This routine visualizes the discrete solutions of the Shallow-Water Equations
%> at the current time step.\n
%> The total height of water can be computed from the 
%> difference of the free surface elevation and the bathymetry. It can be then
%> used to compute approximations of the velocities from the computed momenta
%> via the routine <code>projectDataQ0T2DataDisc()</code>, which requires the 
%> values of the quotients of momenta and height in the quadrature points of 
%> adequate order in each element.\n
%> By this means it is possible to visualize height, free surface elevation, as
%> well as momentum and velocity components for each spatial dimension.\n
%> Furthermore, this routine is responsible for station output, i.e. writing the
%> values of the free surface elevation as well as both velocity components 
%> evaluated at the physical coordinates of each station into a respective
%> field.
%>
%> @param  pd           A struct with problem parameters, precomputed
%>                      fields, and solution data structures, as provided 
%>                      by configureProblem() and preprocessProblem(). 
%>                      @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval pd           The input struct enriched with post-processed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function pd = visualizeSolution(pd, nStep)
varNameElev = {};
varNameVel = {};
vecNames = struct;
dataLagrElev = {};
dataLagrVel = {};

% Evaluate bathymetry (z_b)
if (mod(nStep, pd.outputFrequency(1)) == 0 || pd.isFinished) && any(ismember(pd.outputList, 'bathymetry')) && pd.t >= pd.outputStart(1) && pd.t <= pd.outputEnd(1)
  varNameElev = [ varNameElev, {'bathymetry'} ];
  dataLagrElev = [ dataLagrElev, {projectDataDisc2DataLagr(pd.zbDisc)} ];
end % if

% Evaluate water depth (h)
if (mod(nStep, pd.outputFrequency(1)) == 0 || pd.isFinished) && any(ismember(pd.outputList, 'height')) && pd.t >= pd.outputStart(1) && pd.t <= pd.outputEnd(1)
  varNameElev = [ varNameElev, {'height'} ];
  dataLagrElev = [ dataLagrElev, {projectDataDisc2DataLagr(pd.cDisc(:,:,1))} ];
end % if

% Evaluate water elevation (xi)
if ( (mod(nStep, pd.outputFrequency(1)) == 0 || pd.isFinished) && any(ismember(pd.outputList, 'elevation')) && pd.t >= pd.outputStart(1) && pd.t <= pd.outputEnd(1) ) || ( (mod(nStep, pd.outputFrequency(2)) == 0 || pd.isFinished) && isfield(pd, 'stationElev') && pd.t >= pd.outputStart(2) && pd.t <= pd.outputEnd(2) )
  varNameElev = [ varNameElev, {'elevation'} ];
  dataLagrElev = [ dataLagrElev, {pd.rampInput(pd.t) * projectDataDisc2DataLagr(pd.xiDisc)} ];
end % if

% Save elevation station values
if (mod(nStep, pd.outputFrequency(2)) == 0 || pd.isFinished) && isfield(pd, 'stationElev') && pd.t >= pd.outputStart(2) && pd.t <= pd.outputEnd(2)
  for n = 1 : length(pd.stationElev)
    if pd.p == 0
      dataStationV0T = repmat(dataLagrElev{end}(pd.stationElev{n}(:,1)), 1, 3); % Extract values in vertices of relevant triangles
    else
      dataStationV0T = dataLagrElev{end}(pd.stationElev{n}(:,1),1:3); % Extract values in vertices of relevant triangles
    end % if
    pd.dataElev{n} = [ pd.dataElev{n} ; ...     % Append mean of barycentric weighted values
                       mean(sum(pd.stationElev{n}(:,2:4) .* dataStationV0T, 2)) ];
  end % for
end % if

%% Momentum 
if (mod(nStep, pd.outputFrequency(3)) == 0 || pd.isFinished) && any(ismember(pd.outputList, {'uH', 'vH', 'momentum'})) && pd.t >= pd.outputStart(3) && pd.t <= pd.outputEnd(3)
  vecNames.momentum = {'uH', 'vH'};
  varNameVel = [ varNameVel, vecNames.momentum ];    
  dataLagrVel = [ dataLagrVel, {projectDataDisc2DataLagr(pd.cDisc(:,:,2)), projectDataDisc2DataLagr(pd.cDisc(:,:,3))} ];
end % if

%% Velocity    
% Evaluate velocity
if ( (mod(nStep, pd.outputFrequency(3)) == 0 || pd.isFinished) && any(ismember(pd.outputList, {'u', 'v', 'velocity'}))  && pd.t >= pd.outputStart(3) && pd.t <= pd.outputEnd(3) ) || ( (mod(nStep, pd.outputFrequency(4)) == 0 || pd.isFinished) && isfield(pd, 'stationVel') && pd.t >= pd.outputStart(4) && pd.t <= pd.outputEnd(4) )
  qOrd = max(2*pd.p,1);

  vecNames.velocity = {'u', 'v'};
  varNameVel = [ varNameVel, vecNames.velocity ];
  dataQ0T = (pd.cDisc(:,:,2) * pd.basesOnQuad.phi2D{qOrd}.') ./ (pd.cDisc(:,:,1) * pd.basesOnQuad.phi2D{qOrd}.');
  dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
  dataLagrVel = [ dataLagrVel, {projectDataDisc2DataLagr(dataDisc)} ];

  dataQ0T = (pd.cDisc(:,:,3) * pd.basesOnQuad.phi2D{qOrd}.') ./ (pd.cDisc(:,:,1) * pd.basesOnQuad.phi2D{qOrd}.');
  dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
  dataLagrVel = [ dataLagrVel, {projectDataDisc2DataLagr(dataDisc)} ];
end % if

% Save velocity station values
if (mod(nStep, pd.outputFrequency(4)) == 0 || pd.isFinished) && isfield(pd, 'stationVel') && pd.t >= pd.outputStart(4) && pd.t <= pd.outputEnd(4)
  for n = 1 : length(pd.stationVel)
    if pd.p == 0
      dataStationV0T = repmat(dataLagrVel{end-1}(pd.stationVel{n}(:,1)), 1, 3); % Extract values in vertices of relevant triangles
    else
      dataStationV0T = dataLagrVel{end-1}(pd.stationVel{n}(:,1),1:3); % Extract values in vertices of relevant triangles
    end % if
    pd.dataVel{n,1} = [ pd.dataVel{n,1} ; ...     % Append mean of barycentric weighted values
                        mean(sum(pd.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
    if pd.p == 0
      dataStationV0T = repmat(dataLagrVel{end}(pd.stationVel{n}(:,1)), 1, 3); % Extract values in vertices of relevant triangles
    else
      dataStationV0T = dataLagrVel{end}(pd.stationVel{n}(:,1),1:3); % Extract values in vertices of relevant triangles
    end % if
    pd.dataVel{n,2} = [ pd.dataVel{n,2} ; ...     % Append mean of barycentric weighted values
                        mean(sum(pd.stationVel{n}(:,2:4) .* dataStationV0T, 2)) ];
  end % for
end % if

%% Write visualization output
if pd.outputFrequency(1) == pd.outputFrequency(3) && (~isequal(dataLagrElev, {}) || ~isequal(dataLagrVel, {})) && (mod(nStep, pd.outputFrequency(1)) == 0 || pd.isFinished)
  visualizeDataLagr(pd.g, [dataLagrElev, dataLagrVel], [varNameElev, varNameVel], ['output' filesep pd.name], ceil(nStep / pd.outputFrequency(1)), pd.outputTypes, vecNames);
  if nStep~=0
    pd.changeL2
  end
else
  if ~isequal(dataLagrElev, {}) && (mod(nStep, pd.outputFrequency(1)) || pd.isFinished) == 0
    visualizeDataLagr(pd.g, dataLagrElev, varNameElev, ['output' filesep pd.name '_elev'], ceil(nStep / pd.outputFrequency(1)), pd.outputTypes, struct);
  end % if
  if ~isequal(dataLagrVel, {}) && (mod(nStep, pd.outputFrequency(3))|| pd.isFinished) == 0
    visualizeDataLagr(pd.g, dataLagrVel, varNameVel, ['output' filesep pd.name '_vel'], ceil(nStep / pd.outputFrequency(3)), pd.outputTypes, struct);
  end % if
end % if
end % function
