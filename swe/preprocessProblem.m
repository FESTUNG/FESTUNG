% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file template/preprocessProblem.m
%>
%> @brief Performs all pre-processing tasks, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%===============================================================================
%>
%> @brief Performs all pre-processing steps, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%>
%> This routine is called after template/configureProblem.m.
%>
%> Typically, this step consists of grid generation, computation of derived
%> data structures, pre-computing often needed values (e.g., basis
%> functions on quadrature points), or assembly of time-independent matrix
%> blocks.
%>
%> @param  problemData  A struct with problem parameters, as provided by
%>                      configureProblem(). @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description and precomputed fields.
%>                      @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = preprocessProblem(problemData)
p = problemData.p;

%% Output preprocessing
if ~isdir('output')
  mkdir('output');
end % if

%% Triangulation
switch problemData.gridSource
  case 'debug'
    problemData.g = domainSquare(1);
    problemData.g.idE = zeros(problemData.g.numE,1);
    problemData.g.idE(problemData.g.baryE(:, 2) == 0) = 1; % south
    problemData.g.idE(problemData.g.baryE(:, 1) == 1) = 1; % east
    problemData.g.idE(problemData.g.baryE(:, 2) == 1) = 1; % north
    problemData.g.idE(problemData.g.baryE(:, 1) == 0) = 1; % west
    problemData.g.idE0T = problemData.g.idE(problemData.g.E0T);
    
  case 'hierarchical'
    X1 = [0 1 1 0]; X2 = [0 0 1 1];
    problemData.g = domainHierarchy(X1, X2, problemData.hmax, problemData.refinement);
    problemData.g.idE = zeros(problemData.g.numE,1);
    problemData.g.idE(problemData.g.baryE(:, 2) == 0) = 1; % south
    problemData.g.idE(problemData.g.baryE(:, 1) == 1) = 4; % east
    problemData.g.idE(problemData.g.baryE(:, 2) == 1) = 1; % north
    problemData.g.idE(problemData.g.baryE(:, 1) == 0) = 1; % west
    problemData.g.idE0T = problemData.g.idE(problemData.g.E0T);
    
  case 'ADCIRC'
    projCenter = [problemData.configADCIRC.SLAM0, problemData.configADCIRC.SFEA0];
    [ problemData.g, depth, forcingOS, flowRateRiv ] = domainADCIRC(['swe/fort_' problemData.name '.14'], ...
                                                         ['swe/fort_' problemData.name '.17'], ...
                                                         problemData.configADCIRC.NBFR, ...
                                                         problemData.isSpherical, projCenter);
    
    % Bathymetry
    problemData.zbCont = @(x1,x2) evaluateFuncFromVertexValues(problemData.g, -depth, x1,x2);

    % Convert coordinates to longitude/latitude
    coordSph = [ problemData.g.coordV(:,1) / 6378206.4 / cos(projCenter(2)) + projCenter(1), ...
                 problemData.g.coordV(:,2) / 6378206.4 ];

    % Spatial variation of coriolis parameter
    if problemData.configADCIRC.NCOR == 1 
      problemData.fcCont = @(x1,x2) evaluateFuncFromVertexValues(problemData.g, 2.0 * 7.29212e-5 * sin(coordSph(:,2)), x1,x2);
    else
      problemData.fcCont = @(x1,x2) problemData.configADCIRC.CORI * ones(size(x1));
    end % if
    
    % Newtonian tidal potential
    if problemData.isTidalDomain
      numFrequency = problemData.configADCIRC.NTIF;
      problemData.forcingTidal = cell(2,2,numFrequency);
      problemData.forcingFrequency = cell(2,numFrequency);
      for n = 1 : numFrequency
        problemData.forcingFrequency{1,n} = @(t) cos(problemData.configADCIRC.AMIGT(n) * t);
        problemData.forcingFrequency{2,n} = @(t) -sin(problemData.configADCIRC.AMIGT(n) * t);
        semi = round(0.00014 / problemData.configADCIRC.AMIGT(n));
        switch semi
          case 1
            fX1 = problemData.configADCIRC.TPK(n) * problemData.configADCIRC.FFT(n) * ...
                    problemData.configADCIRC.ETRF(n) * cos(coordSph(:,2)).^2 .* ...
                    cos(pi/180 * (problemData.configADCIRC.FACET(n) + 2*coordSph(:,1)));
            fX2 = problemData.configADCIRC.TPK(n) * problemData.configADCIRC.FFT(n) * ...
                    problemData.configADCIRC.ETRF(n) * cos(coordSph(:,2)).^2 .* ...
                    sin(pi/180 * (problemData.configADCIRC.FACET(n) + 2*coordSph(:,1)));
          case 2
            fX1 = problemData.configADCIRC.TPK(n) * problemData.configADCIRC.FFT(n) * ...
                    problemData.configADCIRC.ETRF(n) * sin(2*coordSph(:,2)) .* ...
                    cos(pi/180 * (problemData.configADCIRC.FACET(n) + coordSph(:,1)));
            fX2 = problemData.configADCIRC.TPK(n) * problemData.configADCIRC.FFT(n) * ...
                    problemData.configADCIRC.ETRF(n) * sin(2*coordSph(:,2)) .* ...
                    sin(pi/180 * (problemData.configADCIRC.FACET(n) + coordSph(:,1)));
          otherwise
            error(['Unsupported tidal potential frequency ' num2str(semi) '!'])
        end % switch
        problemData.forcingTidal{1,1,n} = problemData.gConst ./ (2*problemData.g.areaT) .* ...
          ( (fX1(problemData.g.V0T(:,2)) - fX1(problemData.g.V0T(:,1))) .* problemData.g.B(:,2,2) - ...
            (fX1(problemData.g.V0T(:,3)) - fX1(problemData.g.V0T(:,1))) .* problemData.g.B(:,2,1) );
        problemData.forcingTidal{1,2,n} = problemData.gConst ./ (2*problemData.g.areaT) .* ...
          ( (fX2(problemData.g.V0T(:,2)) - fX2(problemData.g.V0T(:,1))) .* problemData.g.B(:,2,2) - ...
            (fX2(problemData.g.V0T(:,3)) - fX2(problemData.g.V0T(:,1))) .* problemData.g.B(:,2,1) );
        problemData.forcingTidal{2,1,n} = problemData.gConst ./ (2*problemData.g.areaT) .* ...
          (-(fX1(problemData.g.V0T(:,2)) - fX1(problemData.g.V0T(:,1))) .* problemData.g.B(:,1,2) + ...
            (fX1(problemData.g.V0T(:,3)) - fX1(problemData.g.V0T(:,1))) .* problemData.g.B(:,1,1) );
        problemData.forcingTidal{2,2,n} = problemData.gConst ./ (2*problemData.g.areaT) .* ...
          (-(fX2(problemData.g.V0T(:,2)) - fX2(problemData.g.V0T(:,1))) .* problemData.g.B(:,1,2) + ...
            (fX2(problemData.g.V0T(:,3)) - fX2(problemData.g.V0T(:,1))) .* problemData.g.B(:,1,1) );
      end % for
    end % if
    
    % Tidal forcing on open sea boundaries
    numFrequency = problemData.configADCIRC.NBFR;
    problemData.xiFreqOS = cell(2,numFrequency);
    problemData.xiAmpOS = cell(2,numFrequency);
    markTbdrOS = problemData.g.T0E(problemData.g.idE == 4, 1);
    for n = 1 : numFrequency
      problemData.xiFreqOS{1,n} = @(t) cos(problemData.configADCIRC.AMIG(n)*t);
      problemData.xiFreqOS{2,n} = @(t) -sin(problemData.configADCIRC.AMIG(n)*t);
      
      problemData.xiAmpOS{1,n} = sparse(problemData.g.numT,1);
      problemData.xiAmpOS{1,n}(markTbdrOS) = problemData.configADCIRC.FF(n) * forcingOS(n,:,1) .* ...
                                             cos( pi/180 * (problemData.configADCIRC.FACE(n) - forcingOS(n,:,2)) );
      
      problemData.xiAmpOS{2,n} = sparse(problemData.g.numT,1);
      problemData.xiAmpOS{2,n}(markTbdrOS) = problemData.configADCIRC.FF(n) * forcingOS(n,:,1) .* ...
                                             sin( pi/180 * (problemData.configADCIRC.FACE(n) - forcingOS(n,:,2)) );
    end % for
    
    % River inflow
    markEbdrRiv = problemData.g.idE == 3;
    markTbdrRiv = problemData.g.T0E(markEbdrRiv, 1);
    problemData.xiRiv = sparse(problemData.g.numT, 1);
    problemData.xiRiv(markTbdrRiv) = flowRateRiv(:,1);
    problemData.uRiv = sparse(problemData.g.numT, 1);
    problemData.uRiv(markTbdrRiv) = flowRateRiv(:,2) .* problemData.g.nuE(markEbdrRiv,1) - flowRateRiv(:,3) .* problemData.g.nuE(markEbdrRiv,2);
    problemData.vRiv = sparse(problemData.g.numT, 1);
    problemData.vRiv(markTbdrRiv) = flowRateRiv(:,2) .* problemData.g.nuE(markEbdrRiv,2) + flowRateRiv(:,3) .* problemData.g.nuE(markEbdrRiv,1);

    % Stations
    if problemData.isVisStations
      if problemData.configADCIRC.NSTAE == 0 && problemData.configADCIRC.NSTAV == 0
        warning('No stations specified! Disabling station output.')
        problemData.isVisStation = false;
      else
        % Find triangle indices for each station
        coordElev = [ problemData.configADCIRC.XEL, problemData.configADCIRC.YEL ];
        problemData.stationElev = coord2triangle(problemData.g, coordElev(:,1), coordElev(:,2));
        problemData.dataElev = cell(size(problemData.stationElev));
        coordVel = [ problemData.configADCIRC.XEV, problemData.configADCIRC.YEV ];
        problemData.stationVel = coord2triangle(problemData.g, coordVel(:,1), coordVel(:,2));
        problemData.dataVel = cell(size(problemData.stationVel,1),2);
      end % if
    end % if
  otherwise
    error('Invalid gridSource given.')
end % switch

%% Globally constant parameters
problemData.K = problemData.g.numT; % number of triangles
problemData.N = nchoosek(p + 2, p); % number of local DOFs
K = problemData.K;
N = problemData.N;

problemData.g.markE0Tint = problemData.g.idE0T == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrL = problemData.g.idE0T == 1; % [K x 3] mark local edges on the open sea boundary
problemData.g.markE0TbdrRA = problemData.g.idE0T == 2; % [K x 3] mark local edges on the open sea boundary
problemData.g.markE0TbdrRI = problemData.g.idE0T == 3; % [K x 3] mark local edges on the open sea boundary
problemData.g.markE0TbdrOS = problemData.g.idE0T == 4; % [K x 3] mark local edges on the open sea boundary
problemData.isRiverBdr = any(problemData.g.markE0TbdrRI(:));
problemData.g = computeDerivedGridData(problemData.g);

problemData.outputFrequency = max(floor(problemData.numSteps / problemData.outputCount), 1);

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K);

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(N, struct);

%% System matrix for correction of min value exceedence.
problemData.sysMinValueCorrection = [ phi(1,0,0) phi(1,1,0) phi(1,0,1) ; ...
                                      phi(2,0,0) phi(2,1,0) phi(2,0,1) ; ...
                                      phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

%% Computation of matrices on the reference triangle.
refElemPhiPhiPhi = integrateRefElemPhiPhiPhi(N, problemData.basesOnQuad);
refElemDphiLinPhiPhi = integrateRefElemDphiLinPhiPhi(N, problemData.basesOnQuad);
refElemDphiPhi = integrateRefElemDphiPhi(N, problemData.basesOnQuad);
problemData.refElemPhiPhi = integrateRefElemPhiPhi(N, problemData.basesOnQuad);
refEdgePhiIntPhiInt = integrateRefEdgePhiIntPhiInt(N, problemData.basesOnQuad);
refEdgePhiIntPhiExt = integrateRefEdgePhiIntPhiExt(N, problemData.basesOnQuad);
refElemDphiPerQuad = integrateRefElemDphiPerQuad(N, problemData.basesOnQuad);
refEdgePhiIntPerQuad = integrateRefEdgePhiIntPerQuad(N, problemData.basesOnQuad);
problemData.refEdgePhiIntPhiIntPerQuad = integrateRefEdgePhiIntPhiIntPerQuad(N, problemData.basesOnQuad);

%% L2 projections of time-independent algebraic coefficients.
fcDisc = projectFuncCont2DataDisc(problemData.g, problemData.fcCont, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
problemData.zbDisc = projectFuncCont2DataDisc(problemData.g, problemData.zbCont, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);

% Linear representation for globG
zbDiscLin = projectFuncCont2DataDisc(problemData.g, problemData.zbCont, 2, eye(3), computeBasesOnQuad(3, struct)); 

% Evaluate zb in each quadrature point
qOrd = max(2*p,1); [Q, ~] = quadRule1D(qOrd);
problemData.zbPerQuad = cell(3,1);
for n = 1 : 3
  [Q1,Q2] = gammaMap(n, Q);
  problemData.zbPerQuad{n} = problemData.zbCont(problemData.g.mapRef2Phy(1,Q1,Q2), problemData.g.mapRef2Phy(2,Q1,Q2));
end % for

% Visualization of coefficients
if any(ismember(problemData.outputList, 'fc'))
  dataLagr = projectDataDisc2DataLagr(fcDisc);
  visualizeDataLagr(problemData.g, dataLagr, 'f_c', ['output/' problemData.name '_f_c'], 0, problemData.outputTypes);
end % if
if any(ismember(problemData.outputList, 'zb'))
  dataLagr = projectDataDisc2DataLagr(zbDiscLin);
  visualizeDataLagr(problemData.g, dataLagr, 'z_b', ['output/' problemData.name '_z_b'], 0, problemData.outputTypes);
end % if

%% Assembly of time-independent global matrices corresponding to linear contributions.
% Element matrices
globD = assembleMatElemPhiPhiFuncDisc(problemData.g, refElemPhiPhiPhi, fcDisc);
globG = assembleMatElemPhiPhiFuncDiscLin(problemData.g, refElemDphiLinPhiPhi, zbDiscLin);
globH = assembleMatElemDphiPhi(problemData.g, refElemDphiPhi);
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.refElemPhiPhi);

% Edge matrices
globQ = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, problemData.g.areaNuE0Tint);
globQOS = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrOS, refEdgePhiIntPhiInt, problemData.g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrRA, refEdgePhiIntPhiInt, problemData.g.areaNuE0TbdrRA);

% Derived system matrices
problemData.sysW = blkdiag(problemData.globM, problemData.globM, problemData.globM);
problemData.linearTerms = [               sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1},  globQ{2} + globQOS{2} + globQRA{2} - globH{2}; ...
                            problemData.gConst * globG{1},                               sparse(K*N,K*N),                                          -globD; ...
                            problemData.gConst * globG{2},                                         globD,                                 sparse(K*N,K*N) ];
                          
%% Assembly of time-independent global matrices corresponding to non-linear contributions.
% Element matrices
problemData.globF = assembleMatElemDphiPerQuad(problemData.g, refElemDphiPerQuad);

% Edge matrices
[problemData.globRdiag, problemData.globRoffdiag] = assembleMatEdgePhiNuPerQuad(problemData.g, problemData.g.markE0Tint, refEdgePhiIntPerQuad);
problemData.globV = assembleMatEdgePhiPerQuad(problemData.g, refEdgePhiIntPerQuad);

% Boundary matrices
problemData.globRL = assembleMatEdgePhiIntNuPerQuad(problemData.g, problemData.g.markE0TbdrL , refEdgePhiIntPerQuad, problemData.g.areaNuE0TbdrL);
problemData.globROS = assembleMatEdgePhiIntNuPerQuad(problemData.g, problemData.g.markE0TbdrOS, refEdgePhiIntPerQuad, problemData.g.areaNuE0TbdrOS);
globB = assembleMatEdgePhiIntNuPerQuad(problemData.g, problemData.g.markE0TbdrRA, refEdgePhiIntPerQuad, problemData.g.areaNuE0TbdrRA);
problemData.globRdiag = cellfun(@plus, problemData.globRdiag, globB, 'UniformOutput', false);
problemData.globB = assembleMatEdgePhiIntNuPerQuad(problemData.g, problemData.g.markE0TbdrRI, refEdgePhiIntPerQuad, problemData.g.areaNuE0TbdrRI);

%% Assembly of bottom-friction terms.
if problemData.isBottomFrictionVarying
  error('Not implemented')
else
  if problemData.isBottomFrictionNonlinear
    refElemPhiPerQuad = integrateRefElemPhiPerQuad(N, problemData.basesOnQuad);
    problemData.globE = assembleMatElemPhiPhi(problemData.g, refElemPhiPerQuad);
  else
    problemData.globE = problemData.globM;
  end % if
  problemData.globE = problemData.bottomFrictionCoef * problemData.globE;
end % if

%% Assembly of rhs terms.
if problemData.isTidalDomain
  for n = 1 : size(problemData.forcingTidal, 3)
    for i = 1 : 2
      for j = 1 : 2
        problemData.forcingTidal{i,j,n} = assembleMatElemPhiPhiFuncDiscConst(problemData.g, problemData.refElemPhiPhi, problemData.forcingTidal{i,j,n});
      end % for
    end % for
  end % for
end % if
end % function
