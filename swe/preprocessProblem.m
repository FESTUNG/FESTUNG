% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file swe/preprocessProblem.m
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
function pd = preprocessProblem(pd)
%% Triangulation
switch pd.gridSource
  case 'square'
    pd.g = domainSquare(pd.hmax);
    
    % Set edge types
    pd.g.idE = zeros(pd.g.numE,1);
    pd.g.idE(pd.g.baryE(:, 2) == 0) = 3; % south
    pd.g.idE(pd.g.baryE(:, 1) == 1) = 3; % east
    pd.g.idE(pd.g.baryE(:, 2) == 1) = 3; % north
    pd.g.idE(pd.g.baryE(:, 1) == 0) = 3; % west
    
    pd.g.idE0T = pd.g.idE(pd.g.E0T);
    
    % Store edge counts
    pd.g.numEint = sum(pd.g.idE == 0);
    pd.g.numEbdrL = sum(pd.g.idE == 1);
    pd.g.numEbdrRA = sum(pd.g.idE == 2);
    pd.g.numEbdrRI = sum(pd.g.idE == 3);
    pd.g.numEbdrOS = sum(pd.g.idE == 4);
    
    if pd.g.numEbdrOS ~= 0 && pd.isOSCont == 0
      warning('Open sea boundary given but no boundary forcings. Program will use zero boundary condition.')
    end % if
    
    if pd.g.numEbdrRI ~= 0 && pd.isRivCont == 0
      warning('River boundary given but no boundary forcings. Program will use zero boundary condition.')
    end % if
    
  case 'hierarchical'
    X1 = [0 100 100 0]; X2 = [0 0 100 100];
    pd.g = execin('swe/domainHierarchy', X1, X2, pd.hmax, pd.refinement);
    
    % Set edge types
    pd.g.idE = zeros(pd.g.numE,1);
    pd.g.idE(pd.g.baryE(:, 2) == min(X2)) = 3; % south
    pd.g.idE(pd.g.baryE(:, 1) == max(X1)) = 3; % east
    pd.g.idE(pd.g.baryE(:, 2) == max(X2)) = 3; % north
    pd.g.idE(pd.g.baryE(:, 1) == min(X1)) = 3; % west
    pd.g.idE0T = pd.g.idE(pd.g.E0T);
    
    % Store edge counts
    pd.g.numEint = sum(pd.g.idE == 0);
    pd.g.numEbdrL = sum(pd.g.idE == 1);
    pd.g.numEbdrRA = sum(pd.g.idE == 2);
    pd.g.numEbdrRI = sum(pd.g.idE == 3);
    pd.g.numEbdrOS = sum(pd.g.idE == 4);
    
    if pd.g.numEbdrOS ~= 0 && pd.isOSCont == 0
      warning('Open sea boundary given but no boundary forcings. Program will use zero boundary condition.')
    end % if
    
    if pd.g.numEbdrRI ~= 0 && pd.isRivCont == 0
      warning('River boundary given but no boundary forcings. Program will use zero boundary condition.')
    end % if
    
  case 'ADCIRC'
    projCenter = [pd.configADCIRC.SLAM0, pd.configADCIRC.SFEA0];
    [ pd.g, depth, forcingOS, flowRateRiv ] = execin('swe/domainADCIRC', ['swe/fort_' pd.name '.14'], ['swe/fort_' pd.name '.17'], ...
                                                                          pd.configADCIRC.NBFR, pd.isSpherical, projCenter);
    
    % Bathymetry
    h = getFunctionHandle('swe/evaluateFuncFromVertexValues');
    pd.zbCont = @(x1,x2) h(pd.g, -depth, x1,x2);
    
    assert(max( max( abs(depth(pd.g.V0T) + pd.zbCont(pd.g.coordV0T(:,:,1), pd.g.coordV0T(:,:,2))) ) ) < 1.e-5, ...
           'Bathymetry incorrectly constructed!');

    % Convert coordinates to longitude/latitude
    coordSph = [ pd.g.coordV(:,1) / 6378206.4 / cos(projCenter(2)) + projCenter(1), ...
                 pd.g.coordV(:,2) / 6378206.4 ];

    % Spatial variation of coriolis parameter
    if pd.configADCIRC.NCOR == 1 
      pd.fcCont = @(x1,x2) h(pd.g, 2.0 * 7.29212e-5 * sin(coordSph(:,2)), x1,x2);
    else
      pd.fcCont = @(x1,x2) pd.configADCIRC.CORI * ones(size(x1));
    end % if
    
    % Newtonian tidal potential
    if pd.isTidalDomain
      numFrequency = pd.configADCIRC.NTIF;
      pd.forcingTidal = cell(2,2,numFrequency);
      pd.forcingFrequency = cell(2,numFrequency);
      for n = 1 : numFrequency
        pd.forcingFrequency{1,n} = @(t) cos(pd.configADCIRC.AMIGT(n) * t);
        pd.forcingFrequency{2,n} = @(t) -sin(pd.configADCIRC.AMIGT(n) * t);
        semi = round(0.00014 / pd.configADCIRC.AMIGT(n));
        switch semi
          case 1
            fX1 = pd.configADCIRC.TPK(n) * pd.configADCIRC.FFT(n) * ...
                    pd.configADCIRC.ETRF(n) * cos(coordSph(:,2)).^2 .* ...
                    cos(pi/180 * (pd.configADCIRC.FACET(n) + 2*coordSph(:,1)));
            fX2 = pd.configADCIRC.TPK(n) * pd.configADCIRC.FFT(n) * ...
                    pd.configADCIRC.ETRF(n) * cos(coordSph(:,2)).^2 .* ...
                    sin(pi/180 * (pd.configADCIRC.FACET(n) + 2*coordSph(:,1)));
          case 2
            fX1 = pd.configADCIRC.TPK(n) * pd.configADCIRC.FFT(n) * ...
                    pd.configADCIRC.ETRF(n) * sin(2*coordSph(:,2)) .* ...
                    cos(pi/180 * (pd.configADCIRC.FACET(n) + coordSph(:,1)));
            fX2 = pd.configADCIRC.TPK(n) * pd.configADCIRC.FFT(n) * ...
                    pd.configADCIRC.ETRF(n) * sin(2*coordSph(:,2)) .* ...
                    sin(pi/180 * (pd.configADCIRC.FACET(n) + coordSph(:,1)));
          otherwise
            error(['Unsupported tidal potential frequency ' num2str(semi) '!'])
        end % switch
        pd.forcingTidal{1,1,n} = pd.gConst ./ (2*pd.g.areaT) .* ...
          ( (fX1(pd.g.V0T(:,2)) - fX1(pd.g.V0T(:,1))) .* pd.g.B(:,2,2) - ...
            (fX1(pd.g.V0T(:,3)) - fX1(pd.g.V0T(:,1))) .* pd.g.B(:,2,1) );
        pd.forcingTidal{1,2,n} = pd.gConst ./ (2*pd.g.areaT) .* ...
          ( (fX2(pd.g.V0T(:,2)) - fX2(pd.g.V0T(:,1))) .* pd.g.B(:,2,2) - ...
            (fX2(pd.g.V0T(:,3)) - fX2(pd.g.V0T(:,1))) .* pd.g.B(:,2,1) );
        pd.forcingTidal{2,1,n} = pd.gConst ./ (2*pd.g.areaT) .* ...
          (-(fX1(pd.g.V0T(:,2)) - fX1(pd.g.V0T(:,1))) .* pd.g.B(:,1,2) + ...
            (fX1(pd.g.V0T(:,3)) - fX1(pd.g.V0T(:,1))) .* pd.g.B(:,1,1) );
        pd.forcingTidal{2,2,n} = pd.gConst ./ (2*pd.g.areaT) .* ...
          (-(fX2(pd.g.V0T(:,2)) - fX2(pd.g.V0T(:,1))) .* pd.g.B(:,1,2) + ...
            (fX2(pd.g.V0T(:,3)) - fX2(pd.g.V0T(:,1))) .* pd.g.B(:,1,1) );
      end % for
    end % if
    
    % Tidal forcing on open sea boundaries
    if ~pd.isOSCont
      numFrequency = pd.configADCIRC.NBFR;
      pd.xiFreqOS = cell(2,numFrequency);
      pd.xiAmpOS = cell(2,numFrequency);
      markEbdrOS = pd.g.idE == 4;
      for n = 1 : numFrequency
        pd.xiFreqOS{1,n} = @(t) cos(pd.configADCIRC.AMIG(n)*t);
        pd.xiFreqOS{2,n} = @(t) -sin(pd.configADCIRC.AMIG(n)*t);

        xiAmp = zeros(pd.g.numE,1);
        xiAmp(markEbdrOS) = pd.configADCIRC.FF(n) * forcingOS(n,:,1) .* ...
                                               cos( pi/180 * (pd.configADCIRC.FACE(n) - forcingOS(n,:,2)) );
        pd.xiAmpOS{1,n} = xiAmp(pd.g.E0T);

        xiAmp = zeros(pd.g.numE,1);
        xiAmp(markEbdrOS) = pd.configADCIRC.FF(n) * forcingOS(n,:,1) .* ...
                                               sin( pi/180 * (pd.configADCIRC.FACE(n) - forcingOS(n,:,2)) );
        pd.xiAmpOS{2,n} = xiAmp(pd.g.E0T);
      end % for

      if pd.g.numEbdrOS ~= 0 && numFrequency == 0
        warning('Open sea boundary given but no boundary forcings. Program will use zero boundary condition.')
      end % if
    end % if
    
    % River inflow
    if ~pd.isRivCont
      markEbdrRiv = pd.g.idE == 3;
      xiRivE = sparse(pd.g.numE, 1);
      xiRivE(markEbdrRiv) = flowRateRiv(:,1);
      uRivE = sparse(pd.g.numE, 1);
      uRivE(markEbdrRiv) = flowRateRiv(:,2) .* pd.g.nuE(markEbdrRiv,1) - flowRateRiv(:,3) .* pd.g.nuE(markEbdrRiv,2);
      vRivE = sparse(pd.g.numE, 1);
      vRivE(markEbdrRiv) = flowRateRiv(:,2) .* pd.g.nuE(markEbdrRiv,2) + flowRateRiv(:,3) .* pd.g.nuE(markEbdrRiv,1);
      pd.xiRivQ0E0T = xiRivE(pd.g.E0T);
      pd.uRivQ0E0T = uRivE(pd.g.E0T);
      pd.vRivQ0E0T = vRivE(pd.g.E0T);
    end % if

    % Stations
    if pd.isVisStations
      if pd.isAdaptiveTimestep
        error('Station output not implemented for adaptive time stepping.');
      end % if
      if pd.configADCIRC.NSTAE == 0 && pd.configADCIRC.NSTAV == 0
        warning('No stations specified! Disabling station output.')
        pd.isVisStation = false;
      else
        if pd.p > 1
          warning('Note that for station output superlinear solutions are approximated by linear ones.');
        end % if
        % Find triangle indices for each station
        if pd.configADCIRC.NSTAE > 0
          coordElev = [ pd.configADCIRC.XEL, pd.configADCIRC.YEL ];
          pd.stationElev = coord2triangle(pd.g, coordElev(:,1), coordElev(:,2));
          pd.dataElev = cell(size(pd.stationElev));
        end % if
        if pd.configADCIRC.NSTAV > 0
          coordVel = [ pd.configADCIRC.XEV, pd.configADCIRC.YEV ];
          pd.stationVel = coord2triangle(pd.g, coordVel(:,1), coordVel(:,2));
          pd.dataVel = cell(size(pd.stationVel,1),2);
        end % if
      end % if
    end % if
    
    % Clean out ADCIRC config struct
    pd = rmfield(pd, 'configADCIRC');
    
  otherwise
    error('Invalid gridSource given.')
end % switch

if pd.isVisGrid,  visualizeGrid(pd.g);  end % if
%% Globally constant parameters
pd = setdefault(pd, 'outputStart', pd.t0 * ones(1,4));
pd = setdefault(pd, 'outputEnd', pd.tEnd * ones(1,4));
pd = setdefault(pd, 'outputFrequency', max(floor(pd.numSteps / pd.outputCount), 1) * ones(1,4));
pd.K = pd.g.numT; % number of triangles
pd.N = nchoosek(pd.p + 2, pd.p); % number of local DOFs
K = pd.K;
N = pd.N;

pd.g.markE0Tint = pd.g.idE0T == 0; % [K x 3] mark local edges that are interior
pd.g.markE0TbdrL = pd.g.idE0T == 1; % [K x 3] mark local edges on the land boundary
pd.g.markE0TbdrRA = pd.g.idE0T == 2; % [K x 3] mark local edges on the radiation boundary
pd.g.markE0TbdrRI = pd.g.idE0T == 3; % [K x 3] mark local edges on the river boundary
pd.g.markE0TbdrOS = pd.g.idE0T == 4; % [K x 3] mark local edges on the open sea boundary

if ~isempty(pd.slopeLimList)
  pd.g.markV0TbdrRI = ismember(pd.g.V0T, pd.g.V0E(pd.g.E0T(pd.g.markE0TbdrRI), :));
  pd.g.markV0TbdrOS = ismember(pd.g.V0T, pd.g.V0E(pd.g.E0T(pd.g.markE0TbdrOS), :));
  pd.g.markV0TbdrD = pd.g.markV0TbdrRI | pd.g.markV0TbdrOS;
end % if

pd.g = execin('swe/computeDerivedGridData', pd.g);

qOrd1D = 2*pd.p+1; [~, W] = quadRule1D(qOrd1D); numQuad1D = length(W);
pd.g.nuQ0E0T = cell(3,2);
pd.g.nuE0Tprod = cell(3,1);
pd.g.nuE0TsqrDiff = cell(3,1);
pd.g.nuE0Tsqr = cell(3,2);
switch pd.typeBdrL
  case 'natural'
    for n = 1 : 3
      for m = 1 : 2
        pd.g.nuQ0E0T{n,m} = kron(pd.g.nuE0T(:,n,m), ones(numQuad1D, 1));
      end % for
    end % for
  case 'reflected'
    for n = 1 : 3
      pd.g.nuE0Tprod{n} = kron(pd.g.nuE0T(:,n,1) .* pd.g.nuE0T(:,n,2), ones(numQuad1D,1));
      for m = 1 : 2
        pd.g.nuE0Tsqr{n,m} = kron(pd.g.nuE0T(:,n,m) .* pd.g.nuE0T(:,n,m), ones(numQuad1D,1));
        pd.g.nuQ0E0T{n,m} = kron(pd.g.nuE0T(:,n,m), ones(numQuad1D, 1));
      end % for
    end % for
  case 'riemann'
    for n = 1 : 3
      pd.g.nuE0Tprod{n} = kron(pd.g.nuE0T(:,n,1) .* pd.g.nuE0T(:,n,2), ones(numQuad1D,1));
      for m = 1 : 2
        pd.g.nuQ0E0T{n,m} = kron(pd.g.nuE0T(:,n,m), ones(numQuad1D, 1));
      end % for
      pd.g.nuE0TsqrDiff{n} = kron(pd.g.nuE0T(:,n,2) .* pd.g.nuE0T(:,n,2) -pd.g.nuE0T(:,n,1) .* pd.g.nuE0T(:,n,1), ones(numQuad1D,1));
    end % for
  otherwise
    error('Invalid type for land boundary treatment.')
end % switch

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', pd.p, N, K);

%% Lookup table for basis function.
requiredOrders = unique([max(2*pd.p, 1), 2*pd.p+1, 2, 3], 'sorted');
pd.basesOnQuad = computeBasesOnQuad(N, struct, requiredOrders);
if ~isempty(pd.slopeLimList)
  pd.basesOnQuad = computeTaylorBasesV0T(pd.g, N, pd.basesOnQuad);
end % if
basesOnQuadLin = computeBasesOnQuad(3, struct);

%% System matrix for correction of min value exceedence.
pd.sysMinValueCorrection = [ phi(1,0,0) phi(1,1,0) phi(1,0,1) ; ...
                             phi(2,0,0) phi(2,1,0) phi(2,0,1) ; ...
                             phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

%% Computation of matrices on the reference triangle.
pd.refElemPhiPhi = eye(N); %integrateRefElemPhiPhi(N, pd.basesOnQuad);
refElemPhiPhiPhi = execin('swe/integrateRefElemPhiPhiPhi',N, pd.basesOnQuad);
refElemDphiPhi = integrateRefElemDphiPhi(N, pd.basesOnQuad);

refElemPhiLinPhiLin = integrateRefElemPhiPhi(3, basesOnQuadLin);

if pd.p == 0
  refElemPhiPhiDphiLin = permute(integrateRefElemDphiPhiPhi([3 N N], basesOnQuadLin), [3 2 1 4]);
  refElemDphiPhiPhiLin = integrateRefElemDphiPhiPhi([N N 3], basesOnQuadLin);
  refElemPhiPhiPhiLin  = execin('swe/integrateRefElemPhiPhiPhi',[N N 3], basesOnQuadLin);

  refEdgePhiIntPhiIntPhiLin = integrateRefEdgePhiIntPhiIntPhiInt([N N 3], basesOnQuadLin);
  refEdgePhiIntPhiExtPhiLin = permute(execin('swe/integrateRefEdgePhiIntPhiIntPhiExt',[N 3 N], basesOnQuadLin), [1 3 2 4 5]);
else
  refElemPhiPhiDphiLin = permute(integrateRefElemDphiPhiPhi([3 N N], pd.basesOnQuad), [3 2 1 4]);
  refElemDphiPhiPhiLin = integrateRefElemDphiPhiPhi([N N 3], pd.basesOnQuad);
  refElemPhiPhiPhiLin  = execin('swe/integrateRefElemPhiPhiPhi',[N N 3], pd.basesOnQuad);

  refEdgePhiIntPhiIntPhiLin = integrateRefEdgePhiIntPhiIntPhiInt([N N 3], pd.basesOnQuad);
  refEdgePhiIntPhiExtPhiLin = permute(execin('swe/integrateRefEdgePhiIntPhiIntPhiExt',[N 3 N], pd.basesOnQuad), [1 3 2 4 5]);
end % if

refEdgePhiIntPhiInt = integrateRefEdgePhiIntPhiInt(N, pd.basesOnQuad);
refEdgePhiIntPhiExt = integrateRefEdgePhiIntPhiExt(N, pd.basesOnQuad);

refElemDphiPerQuad = execin('swe/integrateRefElemDphiPerQuad', N, pd.basesOnQuad);
refEdgePhiIntPerQuad = execin('swe/integrateRefEdgePhiIntPerQuad', N, pd.basesOnQuad);

%% L2 projections of time-independent algebraic coefficients.
fcDisc = projectFuncCont2DataDisc(pd.g, pd.fcCont, 3, refElemPhiLinPhiLin, basesOnQuadLin);
pd.zbDiscLin = projectFuncCont2DataDisc(pd.g, pd.zbCont, 3, refElemPhiLinPhiLin, basesOnQuadLin);
pd.zbDisc = projectFuncCont2DataDisc(pd.g, pd.zbCont, 2*pd.p+1, pd.refElemPhiPhi, pd.basesOnQuad);
pd.zbLagr = projectDataDisc2DataLagr(pd.zbDiscLin);

% Evaluate zb in each element's quadrature point
[Q1, Q2, ~] = quadRule2D(max(2*pd.p,1)); numQuad2D = length(Q1);
pd.zbQ0T = reshape(pd.zbCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2)).', K * numQuad2D, 1);

% Evaluate zb in each edge's quadrature point
pd.zbQ0E0Tint = cell(3,1);
pd.zbQ0E0Text = cell(3,3);
pd.zbQ0E0TE0T = cell(3,3);
[Q, ~] = quadRule1D(2*pd.p+1); numQuad1D = length(Q);
for nn = 1 : 3
  [Q1, Q2] = gammaMap(nn, Q);
  zbTheta = pd.zbCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2)).';
  pd.zbQ0E0Tint{nn} = reshape(zbTheta, K * numQuad1D, 1);
  for np = 1 : 3
    [QP1, QP2] = theta(nn, np, Q1, Q2);
    zbTheta = pd.zbCont(pd.g.mapRef2Phy(1,QP1,QP2), pd.g.mapRef2Phy(2,QP1,QP2)).';
    pd.zbQ0E0Text{nn,np} = reshape(zbTheta, K * numQuad1D, 1);
    pd.zbQ0E0TE0T{nn,np} = reshape(zbTheta * pd.g.markE0TE0T{nn,np}.', K * numQuad1D, 1);
  end % for
end % for

% Visualization of coefficients
varName = {};
dataLagr = {};
if ismember('f_c', pd.outputList)
  varName = [ varName, {'f_c'} ];
  dataLagr = [ dataLagr, {projectDataDisc2DataLagr(fcDisc)} ];
end % if
if ismember('z_b', pd.outputList)
  varName = [ varName, {'z_b'} ];
  dataLagr = [ dataLagr, {pd.zbLagr} ];
end % if
if ~isempty(varName)
  visualizeDataLagr(pd.g, dataLagr, varName, ['output' filesep pd.name '_coef'], 0, pd.outputTypes);
end % if

%% Assembly of time-independent global matrices corresponding to linear contributions.
% Element matrices
globD = execin('swe/assembleMatElemPhiPhiFuncDisc',pd.g, refElemPhiPhiPhiLin, fcDisc);
globG = execin('swe/assembleMatElemPhiPhiDfuncDisc',pd.g, refElemPhiPhiDphiLin, pd.zbDiscLin);
globH = assembleMatElemDphiPhi(pd.g, refElemDphiPhi);
pd.globM = assembleMatElemPhiPhi(pd.g, pd.refElemPhiPhi);
globU = assembleMatElemDphiPhiFuncDisc(pd.g, refElemDphiPhiPhiLin, pd.zbDiscLin);

% Interior edge matrices
switch pd.typeFlux
  case 'Lax-Friedrichs'
    globQ = assembleMatEdgePhiPhiNu(pd.g, pd.g.markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, pd.g.areaNuE0Tint);
    globO = assembleMatEdgePhiPhiFuncDiscIntNu(pd.g, pd.g.markE0Tint, refEdgePhiIntPhiIntPhiLin, refEdgePhiIntPhiExtPhiLin, pd.zbDiscLin, pd.g.areaNuE0Tint);
  case 'Roe'
    error('not implemented')
    
  otherwise
    error('Unknown flux type')
end % switch

% Boundary edge matrices
globQOS = assembleMatEdgePhiIntPhiIntNu(pd.g, pd.g.markE0TbdrOS, refEdgePhiIntPhiInt, pd.g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu(pd.g, pd.g.markE0TbdrRA, refEdgePhiIntPhiInt, pd.g.areaNuE0TbdrRA);
globOL  = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(pd.g, pd.g.markE0TbdrL, refEdgePhiIntPhiIntPhiLin, pd.zbDiscLin, pd.g.areaNuE0TbdrL);
globORA = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(pd.g, pd.g.markE0TbdrRA, refEdgePhiIntPhiIntPhiLin, pd.zbDiscLin, pd.g.areaNuE0TbdrRA);

% Derived system matrices
pd.sysW = blkdiag(pd.globM, pd.globM, pd.globM);
pd.linearTerms = [ sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1}, globQ{2} + globQOS{2} + globQRA{2} - globH{2}; ...
                   pd.gConst * (globG{1} + globU{1} - globO{1} - globOL{1} - globORA{1}), sparse(K*N,K*N), -globD; ...
                   pd.gConst * (globG{2} + globU{2} - globO{2} - globOL{2} - globORA{2}), globD, sparse(K*N,K*N) ];

% Slope limiting matrices
if ~isempty(pd.slopeLimList)
  globMTaylor = assembleMatElemPhiTaylorPhiTaylor(pd.g, N);
  pd.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(pd.g, N, pd.basesOnQuad);
  pd.globMCorr = spdiags(1 ./ diag(globMTaylor), 0, K*N, K*N) * globMTaylor;
end % if

%% Assembly of time-independent global matrices corresponding to non-linear contributions.
% Element matrices
pd.globF = execin('swe/assembleMatElemDphiPerQuad',pd.g, refElemDphiPerQuad);

% Edge matrices
[pd.globRdiag, pd.globRoffdiag] = execin('swe/assembleMatEdgePhiNuPerQuad',pd.g, pd.g.markE0Tint, refEdgePhiIntPerQuad);
pd.globV = execin('swe/assembleMatEdgePhiPerQuad',pd.g, refEdgePhiIntPerQuad);

% Bottom-friction terms
if pd.isBottomFrictionVarying
  bottomFrictionDisc = projectFuncCont2DataDisc(pd.g, pd.bottomFrictionCont, 2*pd.p, pd.refElemPhiPhi);
  if pd.isBottomFrictionNonlinear
    refElemPhiPhiPerQuad = integrateRefElemPhiPhiPerQuad(N, pd.basesOnQuad);
    pd.globE = assembleMatElemPhiFuncDiscPerQuad(pd.g, refElemPhiPhiPerQuad, bottomFrictionDisc);
  else
    pd.globE = execin('swe/assembleMatElemPhiPhiFuncDisc',pd.g, refElemPhiPhiPhi, bottomFrictionDisc);
  end % if
else
  if pd.isBottomFrictionNonlinear
    % Assembling the non linear bottom friction matrix globE, which will then be
    % applied to each component of |u|*u evaluated in each quadrature point.
    % Thus, globE is a block-diagonal matrix of size KxR, assembled from a
    % reference block of size NxR. The assembly operation itself is just
    % the Kronecker product of a diagonal matrix with the transformation
    % determinants 2*|T| and the reference block - precisely the same
    % operation as for the assembly of a mass matrix. The corresponding
    % routine assembleMatElemPhiPhi is therefore re-used here.
    refElemPhiPerQuad = execin('swe/integrateRefElemPhiPerQuad',N, pd.basesOnQuad);
    pd.globE = pd.bottomFrictionCoef * assembleMatElemPhiPhi(pd.g, refElemPhiPerQuad);
  else
    pd.globE = pd.bottomFrictionCoef * pd.globM;
  end % if
end % if

% Boundary matrices
if pd.g.numEbdrL > 0 % Land boundaries
  pd.globRL = execin('swe/assembleMatEdgePhiIntNuPerQuad',pd.g, pd.g.markE0TbdrL, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrL);
  if strcmp(pd.typeBdrL, 'riemann')
      pd.globVL = execin('swe/assembleMatEdgePhiIntPerQuad',pd.g, pd.g.markE0TbdrL, refEdgePhiIntPerQuad, pd.g.areaE0TbdrL);
  end % if
end % if

if pd.g.numEbdrRA > 0 % Radiation boundaries
  globRRA = execin('swe/assembleMatEdgePhiIntNuPerQuad',pd.g, pd.g.markE0TbdrRA, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrRA);
  pd.globRdiag = cellfun(@plus, pd.globRdiag, globRRA, 'UniformOutput', false);
end % if

pd.globLRI = { sparse(K*N,1); sparse(K*N,1); sparse(K*N,1) };
if pd.g.numEbdrRI > 0 % River boundaries
  if ~pd.isRivCont
    if any(ismember(pd.slopeLimList, 'momentum'))
      pd.zbV0T = computeFuncContV0T(pd.g, pd.zbCont);

      % To determine the vertex values for all triangles we first compute
      % the vertex values for triangles that have a boundary edge of river
      % type. We average the values from both sides of the vertex if the 
      % triangle has more than one boundary edge.
      xiV0T = [ sum(pd.xiRivQ0E0T(:,[2,3]),2) ./ sum(pd.g.markE0TbdrRI(:,[2,3]),2), ...
                sum(pd.xiRivQ0E0T(:,[1,3]),2) ./ sum(pd.g.markE0TbdrRI(:,[1,3]),2), ...
                sum(pd.xiRivQ0E0T(:,[1,2]),2) ./ sum(pd.g.markE0TbdrRI(:,[1,2]),2) ];
      uV0T = [ sum(pd.uRivQ0E0T(:,[2,3]),2) ./ sum(pd.g.markV0TbdrRI(:,[2,3]),2), ...
               sum(pd.uRivQ0E0T(:,[1,3]),2) ./ sum(pd.g.markV0TbdrRI(:,[1,3]),2), ...
               sum(pd.uRivQ0E0T(:,[1,2]),2) ./ sum(pd.g.markV0TbdrRI(:,[1,2]),2) ];
      vV0T = [ sum(pd.vRivQ0E0T(:,[2,3]),2) ./ sum(pd.g.markV0TbdrRI(:,[2,3]),2), ...
               sum(pd.vRivQ0E0T(:,[1,3]),2) ./ sum(pd.g.markV0TbdrRI(:,[1,3]),2), ...
               sum(pd.vRivQ0E0T(:,[1,2]),2) ./ sum(pd.g.markV0TbdrRI(:,[1,2]),2) ];

      % all vertex indices for which values have to be set
      nonUniqueVertices = pd.g.V0T(pd.g.markV0TbdrRI);
      % the vertex values for triangles with boundary edges and (possibly) 
      % NaN if it is not a boundary vertex
      dataV0T = [ xiV0T(pd.g.markV0TbdrRI), uV0T(pd.g.markV0TbdrRI), vV0T(pd.g.markV0TbdrRI) ];
      
      % add the values of all contributing vertices via matrix-vector 
      % product
      vertInd2VertIndUniqueRI = bsxfun(@eq, (1:pd.g.numV).', nonUniqueVertices.');
      
      dataV = zeros(pd.g.numV,3);
      dataVCountRI = zeros(pd.g.numV,3);
      
      dataV(:,1) = vertInd2VertIndUniqueRI * setNaN2Zero(dataV0T(:,1));
      dataVCountRI(:,1) = vertInd2VertIndUniqueRI * double(~isnan(dataV0T(:,1)));
      dataV(:,2) = vertInd2VertIndUniqueRI * setNaN2Zero(dataV0T(:,2));
      dataVCountRI(:,2) = vertInd2VertIndUniqueRI * double(~isnan(dataV0T(:,2)));
      dataV(:,3) = vertInd2VertIndUniqueRI * setNaN2Zero(dataV0T(:,3));
      dataVCountRI(:,3) = vertInd2VertIndUniqueRI * double(~isnan(dataV0T(:,3)));
      
      % dataV contains the sum of all triangles vertex values and is
      % divided by the number of contributing elements in dataVCount
      dataV = dataV ./ dataVCountRI;
      
      xiV = dataV(:, 1);
      if any(ismember(pd.slopeLimList, 'elevation'))
        pd.xiV0Triv = xiV(pd.g.V0T);
      end % if
      hV0T = xiV(pd.g.V0T) - pd.zbV0T;
      uV = dataV(:, 2);
      pd.uHV0Triv = uV(pd.g.V0T) .* hV0T;
      vV = dataV(:, 3);
      pd.vHV0Triv = vV(pd.g.V0T) .* hV0T;
    elseif any(ismember(pd.slopeLimList, 'elevation'))
      % To determine the vertex values for all triangles we first compute
      % the vertex values for triangles that have a boundary edge of river
      % type. We average the values from both sides of the vertex if the 
      % triangle has more than one boundary edge.
      xiV0T = [ sum(pd.xiRivQ0E0T(:,[2,3]),2) ./ sum(pd.g.markE0TbdrRI(:,[2,3]),2), ...
                sum(pd.xiRivQ0E0T(:,[1,3]),2) ./ sum(pd.g.markE0TbdrRI(:,[1,3]),2), ...
                sum(pd.xiRivQ0E0T(:,[1,2]),2) ./ sum(pd.g.markE0TbdrRI(:,[1,2]),2) ];
      
      % all vertex indices for which values have to be set
      nonUniqueVertices = pd.g.V0T(pd.g.markV0TbdrRI);
      % the vertex values for triangles with boundary edges and (possibly) 
      % NaN if it is not a boundary vertex
      xiV0T = xiV0T(pd.g.markV0TbdrRI);
      
      % add the values of all contributing vertices via matrix-vector 
      % product
      vertInd2VertIndUniqueRI = bsxfun(@eq, (1:pd.g.numV).', nonUniqueVertices.');
      
      xiV = vertInd2VertIndUniqueRI * setNaN2Zero(xiV0T);
      xiVCountRI = vertInd2VertIndUniqueRI * double(~isnan(xiV0T));
      
      % xiV contains the sum of all triangles vertex values and is divided
      % by the number of contributing elements in dataVCount
      xiV = xiV ./ xiVCountRI;
      pd.xiV0Triv = xiV(pd.g.V0T);
    end % if

    pd.xiRivQ0E0T = kron(pd.xiRivQ0E0T, ones(numQuad1D,1));
    pd.uRivQ0E0T = kron(pd.uRivQ0E0T, ones(numQuad1D,1));
    pd.vRivQ0E0T = kron(pd.vRivQ0E0T, ones(numQuad1D,1));

  elseif any(ismember(pd.slopeLimList, 'momentum'))
    pd.zbV0T = computeFuncContV0T(pd.g, pd.zbCont);
  end % if
  
  pd.globRRI = execin('swe/assembleMatEdgePhiIntNuPerQuad',pd.g, pd.g.markE0TbdrRI, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrRI);
  
  if pd.isRiemRiv
    pd.globVRI = execin('swe/assembleMatEdgePhiIntPerQuad',pd.g, pd.g.markE0TbdrRI, refEdgePhiIntPerQuad, pd.g.areaE0TbdrRI);
  end % if
  
  if ~pd.isRamp && ~pd.isRivCont && ~pd.isRiemRiv
    if pd.isCoupling
      pd.massFluxQ0E0TRiv = zeros(K,3,numQuad1D);
    end % if
    for n = 1 : 3
      hRiv = pd.xiRivQ0E0T(:,n) - pd.zbQ0E0Tint{n};
      uHRiv = pd.uRivQ0E0T(:,n) .* hRiv;
      vHRiv = pd.vRivQ0E0T(:,n) .* hRiv;
      
      uuHRiv = pd.uRivQ0E0T(:,n) .* uHRiv;
      uvHRiv = pd.uRivQ0E0T(:,n) .* vHRiv;
      vvHRiv = pd.vRivQ0E0T(:,n) .* vHRiv;
      gHHRiv = pd.gConst * pd.xiRivQ0E0T(:,n) .* ( 0.5 * pd.xiRivQ0E0T(:,n) - pd.zbQ0E0Tint{n} );
      
      pd.globLRI{1} = pd.globLRI{1} + pd.globRRI{n,1} * uHRiv + pd.globRRI{n,2} * vHRiv;
      pd.globLRI{2} = pd.globLRI{2} + pd.globRRI{n,1} * (uuHRiv + gHHRiv) + pd.globRRI{n,2} * uvHRiv;
      pd.globLRI{3} = pd.globLRI{3} + pd.globRRI{n,1} * uvHRiv + pd.globRRI{n,2} * (vvHRiv + gHHRiv);
      
      if pd.isCoupling
        pd.massFluxQ0E0TRiv(:,n,:) = bsxfun(@times, reshape(uHRiv.*pd.g.nuQ0E0T{n,1}+vHRiv.*pd.g.nuQ0E0T{n,2}, [numQuad1D, K])', pd.g.markE0TbdrRI(:,n));
      end % if
    end % for
  end % if
else
  if any(ismember(pd.slopeLimList, 'elevation'))
    pd.xiV0Triv = sparse(K,3);
  end % if
  if any(ismember(pd.slopeLimList, 'momentum'))
    pd.uHV0Triv = sparse(K,3);
    pd.vHV0Triv = sparse(K,3);
  end % if
end % if

if pd.g.numEbdrOS > 0 % Open sea boundaries
  if ~pd.isOSCont && any(ismember(pd.slopeLimList, 'elevation'))
    xiE0Tos = zeros(K,3);
    for n = 1 : numFrequency
      xiE0Tos = xiE0Tos + pd.xiFreqOS{1,n}(pd.t0) * pd.xiAmpOS{1,n} + pd.xiFreqOS{2,n}(pd.t0) * pd.xiAmpOS{2,n};
    end % for
    % To determine the vertex values for all triangles we first compute
    % the vertex values for triangles that have a boundary edge of open sea
    % type. We average the values from both sides of the vertex if the 
    % triangle has more than one boundary edge.
    xiV0T = [ sum(xiE0Tos(:,[2,3]),2) ./ sum(pd.g.markE0TbdrOS(:,[2,3]),2), ...
              sum(xiE0Tos(:,[1,3]),2) ./ sum(pd.g.markE0TbdrOS(:,[1,3]),2), ...
              sum(xiE0Tos(:,[1,2]),2) ./ sum(pd.g.markE0TbdrOS(:,[1,2]),2) ];

    % all vertex indices for which values have to be set
    nonUniqueVertices = pd.g.V0T(pd.g.markV0TbdrOS);
    % the vertex values for triangles with boundary edges and (possibly) 
    % NaN if it is not a boundary vertex
    xiV0T = xiV0T(pd.g.markV0TbdrOS);
    
    % add the values of all contributing vertices via matrix-vector product
    pd.vertInd2VertIndUniqueOS = bsxfun(@eq, (1:pd.g.numV).', nonUniqueVertices.');
    
    xiV = pd.vertInd2VertIndUniqueOS * setNaN2Zero(xiV0T);
    pd.xiVCountOS = pd.vertInd2VertIndUniqueOS * double(~isnan(xiV0T));
    
    % xiV contains the sum of the vertex values of all triangles and is
    % divided by the number of contributing elements in dataVCountOS
    xiV = xiV ./ pd.xiVCountOS;
    pd.xiV0Tos = xiV(pd.g.V0T);
  end % if
  
  pd.globROS = execin('swe/assembleMatEdgePhiIntNuPerQuad',pd.g, pd.g.markE0TbdrOS, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrOS);
  if pd.isRiemOS
    pd.globVOS = execin('swe/assembleMatEdgePhiIntPerQuad',pd.g, pd.g.markE0TbdrOS, refEdgePhiIntPerQuad, pd.g.areaE0TbdrOS);
  end % if
elseif any(ismember(pd.slopeLimList, 'elevation'))
  pd.xiV0Tos = sparse(K,3);
end % if

if pd.isRamp && (pd.isRivCont || pd.isOSCont)
  error('Ramping is not supported for algebraic boundary conditions.');
end % 

%% Assembly of rhs terms.
% Assemble Newtonian tide potential matrix
if pd.isTidalDomain
  if pd.p == 0
    refElemPhiPhiConstPhiLeastLin = execin('swe/integrateRefElemPhiPhiPhi', [N 3 1], basesOnQuadLin);
  else
    refElemPhiPhiConstPhiLeastLin = execin('swe/integrateRefElemPhiPhiPhi', [N N 1], pd.basesOnQuad);
  end % if
  for n = 1 : size(pd.forcingTidal, 3)
    for i = 1 : 2
      for j = 1 : 2
        pd.forcingTidal{i,j,n} = assembleMatElemPhiPhiFuncDisc(pd.g, refElemPhiPhiConstPhiLeastLin, pd.forcingTidal{i,j,n});
      end % for
    end % for
  end % for
end % if

%% Variable timestepping preparation.
if pd.isAdaptiveTimestep
  if pd.p > 1
    warning('Note that for time step adaptivity superlinear solutions are approximated by linear ones.');
  end % if
  pd.avgDiff = sum(abs(pd.g.coordV0T(:,[1 2 3],:) - pd.g.coordV0T(:,[2 3 1],:)), 2) / 3;
  pd.avgDepth = -sum(pd.zbCont(pd.g.coordV0T(:,:,1), pd.g.coordV0T(:,:,2)), 2) / 3;
end % if

%% Function handles
pd.swe_visualizeSolution = getFunctionHandle('swe/visualizeSolution');
pd.computeAveragedVariablesQ0E0Tint = getFunctionHandle('swe/computeAveragedVariablesQ0E0Tint');
pd.computeAveragedVariablesQ0E0Tland = getFunctionHandle('swe/computeAveragedVariablesQ0E0Tland');
pd.computeAveragedVariablesQ0E0Triv = getFunctionHandle('swe/computeAveragedVariablesQ0E0Triv');
pd.computeAveragedVariablesQ0E0Tos = getFunctionHandle('swe/computeAveragedVariablesQ0E0Tos');
pd.computeLaxFriedrichsCoefficient = getFunctionHandle('swe/computeLaxFriedrichsCoefficient');
end % function
