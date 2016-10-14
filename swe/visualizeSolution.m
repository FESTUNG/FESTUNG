% Auxiliary function for visualization of 2D Shallow-Water Equations output.

%===============================================================================
%> @file visualizeSolution.m
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
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>                      Modified 08/23/16 by Hennes Hajduk
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
if mod(nStep, pd.outputFrequency) == 0 || pd.isFinished
  nOutput = ceil(nStep / pd.outputFrequency);
  varName = {};
  vecNames = struct;
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
    vecNames.momentum = {'uH', 'vH'};
    varName = [ varName, vecNames.momentum ];    
    dataLagr = [ dataLagr, {projectDataDisc2DataLagr(pd.cDisc(:,:,2)), projectDataDisc2DataLagr(pd.cDisc(:,:,3))} ];
  end % if
  
  %% Velocity    
  % Evaluate velocity
  if any(ismember(pd.outputList, {'u', 'v', 'velocity'})) || isfield(pd, 'stationVel')
    vecNames.velocity = {'u', 'v'};
    varName = [ varName, vecNames.velocity ];
    dataQ0T = (pd.cDisc(:,:,2) * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.') ./ (hDisc * pd.basesOnQuad.phi2D{max(2*pd.p,1)}.');
    dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
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
    dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*pd.p, pd.refElemPhiPhi, pd.basesOnQuad);
    dataLagr = projectDataDisc2DataLagr(dataDisc);
  end % if

  % Visualize y-velocity (v)
  if any(ismember(pd.outputList, 'v'))
    visualizeDataLagr(pd.g, dataLagr, 'v_h', ['output/' pd.name '_v'], nOutput, pd.outputTypes);
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
  visualizeDataLagr(pd.g, dataLagr, varName, ['output' filesep pd.name], nOutput, pd.outputTypes, vecNames);
end % if
end % function

