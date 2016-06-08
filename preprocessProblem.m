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
function pd = preprocessProblem(pd)
%% Triangulation
switch pd.gridSource
  case 'square'
    pd.g = domainSquare(pd.hmax);
    
    % Set edge types
    pd.g.idE = zeros(pd.g.numE,1);
    pd.g.idE(pd.g.baryE(:, 2) == 0) = 1; % south
    pd.g.idE(pd.g.baryE(:, 1) == 1) = 4; % east
    pd.g.idE(pd.g.baryE(:, 2) == 1) = 1; % north
    pd.g.idE(pd.g.baryE(:, 1) == 0) = 1; % west
    pd.g.idE0T = pd.g.idE(pd.g.E0T);
    
    % Store edge counts
    pd.g.numEint = sum(pd.g.idE == 0);
    pd.g.numEbdrL = sum(pd.g.idE == 1);
    pd.g.numEbdrRA = sum(pd.g.idE == 2);
    pd.g.numEbdrRI = sum(pd.g.idE == 3);
    pd.g.numEbdrOS = sum(pd.g.idE == 4);
    
  case 'hierarchical'
    X1 = [0 1 1 0]; X2 = [0 0 1 1];
    pd.g = domainHierarchy(X1, X2, pd.hmax, pd.refinement);
    
    % Set edge types
    pd.g.idE = zeros(pd.g.numE,1);
    pd.g.idE(pd.g.baryE(:, 2) == 0) = 1; % south
    pd.g.idE(pd.g.baryE(:, 1) == 1) = 4; % east
    pd.g.idE(pd.g.baryE(:, 2) == 1) = 1; % north
    pd.g.idE(pd.g.baryE(:, 1) == 0) = 1; % west
    pd.g.idE0T = pd.g.idE(pd.g.E0T);
    
    % Store edge counts
    pd.g.numEint = sum(pd.g.idE == 0);
    pd.g.numEbdrL = sum(pd.g.idE == 1);
    pd.g.numEbdrRA = sum(pd.g.idE == 2);
    pd.g.numEbdrRI = sum(pd.g.idE == 3);
    pd.g.numEbdrOS = sum(pd.g.idE == 4);
    
  case 'ADCIRC'
    projCenter = [pd.configADCIRC.SLAM0, pd.configADCIRC.SFEA0];
    [ pd.g, depth, forcingOS, flowRateRiv ] = domainADCIRC(['swe/fort_' pd.name '.14'], ['swe/fort_' pd.name '.17'], ...
                                                           pd.configADCIRC.NBFR, pd.isSpherical, projCenter);
    
    % Bathymetry
    pd.zbCont = @(x1,x2) evaluateFuncFromVertexValues(pd.g, -depth, x1,x2);
    assert(max( max( abs(depth(pd.g.V0T) + pd.zbCont(pd.g.coordV0T(:,:,1), pd.g.coordV0T(:,:,2))) ) ) < 1.e-5, ...
           'Bathymetry incorrectly constructed!');

    % Convert coordinates to longitude/latitude
    coordSph = [ pd.g.coordV(:,1) / 6378206.4 / cos(projCenter(2)) + projCenter(1), ...
                 pd.g.coordV(:,2) / 6378206.4 ];

    % Spatial variation of coriolis parameter
    if pd.configADCIRC.NCOR == 1 
      pd.fcCont = @(x1,x2) evaluateFuncFromVertexValues(pd.g, 2.0 * 7.29212e-5 * sin(coordSph(:,2)), x1,x2);
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
    numFrequency = pd.configADCIRC.NBFR;
    pd.xiFreqOS = cell(2,numFrequency);
    pd.xiAmpOS = cell(2,numFrequency);
    markTbdrOS = pd.g.T0E(pd.g.idE == 4, 1);
    for n = 1 : numFrequency
      pd.xiFreqOS{1,n} = @(t) cos(pd.configADCIRC.AMIG(n)*t);
      pd.xiFreqOS{2,n} = @(t) -sin(pd.configADCIRC.AMIG(n)*t);
      
      pd.xiAmpOS{1,n} = sparse(pd.g.numT,1);
      pd.xiAmpOS{1,n}(markTbdrOS) = pd.configADCIRC.FF(n) * forcingOS(n,:,1) .* ...
                                             cos( pi/180 * (pd.configADCIRC.FACE(n) - forcingOS(n,:,2)) );
      
      pd.xiAmpOS{2,n} = sparse(pd.g.numT,1);
      pd.xiAmpOS{2,n}(markTbdrOS) = pd.configADCIRC.FF(n) * forcingOS(n,:,1) .* ...
                                             sin( pi/180 * (pd.configADCIRC.FACE(n) - forcingOS(n,:,2)) );
    end % for
    
    % River inflow
    markEbdrRiv = pd.g.idE == 3;
    markTbdrRiv = pd.g.T0E(markEbdrRiv, 1);
    pd.xiRiv = sparse(pd.g.numT, 1);
    pd.xiRiv(markTbdrRiv) = flowRateRiv(:,1);
    pd.uRiv = sparse(pd.g.numT, 1);
    pd.uRiv(markTbdrRiv) = flowRateRiv(:,2) .* pd.g.nuE(markEbdrRiv,1) - flowRateRiv(:,3) .* pd.g.nuE(markEbdrRiv,2);
    pd.vRiv = sparse(pd.g.numT, 1);
    pd.vRiv(markTbdrRiv) = flowRateRiv(:,2) .* pd.g.nuE(markEbdrRiv,2) + flowRateRiv(:,3) .* pd.g.nuE(markEbdrRiv,1);

    % Stations
    if pd.isVisStations
      error('not implemented')
      if pd.configADCIRC.NSTAE == 0 && pd.configADCIRC.NSTAV == 0
        warning('No stations specified! Disabling station output.')
        pd.isVisStation = false;
      else
        % Find triangle indices for each station
        coordElev = [ pd.configADCIRC.XEL, pd.configADCIRC.YEL ];
        pd.stationElev = coord2triangle(pd.g, coordElev(:,1), coordElev(:,2));
        pd.dataElev = cell(size(pd.stationElev));
        coordVel = [ pd.configADCIRC.XEV, pd.configADCIRC.YEV ];
        pd.stationVel = coord2triangle(pd.g, coordVel(:,1), coordVel(:,2));
        pd.dataVel = cell(size(pd.stationVel,1),2);
      end % if
    end % if
    
    % Clean out ADCIRC config struct
    pd = rmfield(pd, pd.configADCIRC);
    
  otherwise
    error('Invalid gridSource given.')
end % switch

%% Globally constant parameters
pd.outputFrequency = max(floor(pd.numSteps / pd.outputCount), 1);
pd.K = pd.g.numT; % number of triangles
pd.N = nchoosek(pd.p + 2, pd.p); % number of local DOFs
K = pd.K;
N = pd.N;

pd.g.markE0Tint = pd.g.idE0T == 0; % [K x 3] mark local edges that are interior
pd.g.markE0TbdrL = pd.g.idE0T == 1; % [K x 3] mark local edges on the open sea boundary
pd.g.markE0TbdrRA = pd.g.idE0T == 2; % [K x 3] mark local edges on the open sea boundary
pd.g.markE0TbdrRI = pd.g.idE0T == 3; % [K x 3] mark local edges on the open sea boundary
pd.g.markE0TbdrOS = pd.g.idE0T == 4; % [K x 3] mark local edges on the open sea boundary

pd.g = computeDerivedGridData(pd.g);

qOrd1D = 2*pd.p+1; [~, W] = quadRule1D(qOrd1D); numQuad1D = length(W);
pd.g.nuQ0E0T = cell(3,2);
pd.g.nuE0Tprod = cell(3,1);
pd.g.nuE0TsqrDiff = cell(3,1);
pd.g.nuE0Tsqr = cell(3,2);
for n = 1 : 3
  pd.g.nuE0Tprod{n} = kron(pd.g.nuE0T(:,n,1) .* pd.g.nuE0T(:,n,2), ones(numQuad1D,1));
  for m = 1 : 2
    pd.g.nuE0Tsqr{n,m} = kron(pd.g.nuE0T(:,n,m) .* pd.g.nuE0T(:,n,m), ones(numQuad1D,1));
    pd.g.nuQ0E0T{n,m} = kron(pd.g.nuE0T(:,n,m), ones(numQuad1D, 1));
  end % for
  pd.g.nuE0TsqrDiff{n} = pd.g.nuE0Tsqr{n,2} - pd.g.nuE0Tsqr{n,1};
end % for

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', pd.p, N, K);

%% Lookup table for basis function.
pd.basesOnQuad = computeBasesOnQuad(N, struct);
if ~isempty(pd.slopeLimList)
  pd.basesOnQuad = computeTaylorBasesV0T(pd.g, N, pd.basesOnQuad);
end % if

%% System matrix for correction of min value exceedence.
pd.sysMinValueCorrection = [ phi(1,0,0) phi(1,1,0) phi(1,0,1) ; ...
                             phi(2,0,0) phi(2,1,0) phi(2,0,1) ; ...
                             phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

%% Computation of matrices on the reference triangle.
pd.refElemPhiPhi = integrateRefElemPhiPhi(N, pd.basesOnQuad);
refElemPhiPhiPhi = integrateRefElemPhiPhiPhi(N, pd.basesOnQuad);
refElemDphiPhi = integrateRefElemDphiPhi(N, pd.basesOnQuad);

refElemDphiLinPhiPhi = integrateRefElemDphiPhiPhi([3 N N], pd.basesOnQuad);
refElemDphiPhiPhiLin = integrateRefElemDphiPhiPhi([N N 3], pd.basesOnQuad);
refElemPhiPhiPhiLin  = integrateRefElemPhiPhiPhi([N N 3], pd.basesOnQuad);

refEdgePhiIntPhiInt = integrateRefEdgePhiIntPhiInt(N, pd.basesOnQuad);
refEdgePhiIntPhiExt = integrateRefEdgePhiIntPhiExt(N, pd.basesOnQuad);

refEdgePhiIntPhiIntPhiLin = integrateRefEdgePhiIntPhiIntPhiInt([N N 3], pd.basesOnQuad);
refEdgePhiIntPhiExtPhiLin = permute(integrateRefEdgePhiIntPhiIntPhiExt([N 3 N], pd.basesOnQuad), [1 3 2 4 5]);

refElemDphiPerQuad = integrateRefElemDphiPerQuad(N, pd.basesOnQuad);
refEdgePhiIntPerQuad = integrateRefEdgePhiIntPerQuad(N, pd.basesOnQuad);
pd.refEdgePhiIntPhiIntPerQuad = integrateRefEdgePhiIntPhiIntPerQuad(N, pd.basesOnQuad);

%% L2 projections of time-independent algebraic coefficients.
basesOnQuadLin = computeBasesOnQuad(3, struct);
refElemPhiLinPhiLin = integrateRefElemPhiPhi(3, basesOnQuadLin);

fcDisc = projectFuncCont2DataDisc(pd.g, pd.fcCont, 2*pd.p, refElemPhiLinPhiLin, basesOnQuadLin);
pd.zbDisc = projectFuncCont2DataDisc(pd.g, pd.zbCont, 2*pd.p, refElemPhiLinPhiLin, basesOnQuadLin);

% Linear representation for globG
zbDiscLin = projectFuncCont2DataDisc(pd.g, pd.zbCont, 2, eye(3), basesOnQuadLin); 

% Evaluate zb in each quadrature point
qOrd = max(2*pd.p,1); [Q, ~] = quadRule1D(qOrd);
pd.zbPerQuad = cell(3,1);
for n = 1 : 3
  [Q1,Q2] = gammaMap(n, Q);
  pd.zbPerQuad{n} = pd.zbCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2));
end % for

if ~isempty(pd.slopeLimList)
  pd.zbV0T = computeFuncContV0T(pd.g, pd.zbCont);
end % if

% Visualization of coefficients
if any(ismember(pd.outputList, 'fc'))
  dataLagr = projectDataDisc2DataLagr(fcDisc);
  visualizeDataLagr(pd.g, dataLagr, 'f_c', ['output/' pd.name '_f_c'], 0, pd.outputTypes);
end % if
if any(ismember(pd.outputList, 'zb'))
  dataLagr = projectDataDisc2DataLagr(zbDiscLin);
  visualizeDataLagr(pd.g, dataLagr, 'z_b', ['output/' pd.name '_z_b'], 0, pd.outputTypes);
end % if

%% Assembly of time-independent global matrices corresponding to linear contributions.
% Element matrices
globD = assembleMatElemPhiPhiFuncDisc(pd.g, refElemPhiPhiPhiLin, fcDisc);
globG = assembleMatElemPhiPhiDfuncDisc(pd.g, refElemDphiLinPhiPhi, zbDiscLin);
globH = assembleMatElemDphiPhi(pd.g, refElemDphiPhi);
pd.globM = assembleMatElemPhiPhi(pd.g, pd.refElemPhiPhi);
globU = assembleMatElemDphiPhiFuncDisc(pd.g, refElemDphiPhiPhiLin, zbDiscLin);

% Interior edge matrices
switch pd.typeFlux
  case 'Roe'
    error('not implemented')
    
  case 'Lax-Friedrichs'
    globQ = assembleMatEdgePhiPhiNu(pd.g, pd.g.markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, pd.g.areaNuE0Tint);
    globO = assembleMatEdgePhiPhiFuncDiscNu(pd.g, pd.g.markE0Tint, refEdgePhiIntPhiIntPhiLin, refEdgePhiIntPhiExtPhiLin, zbDiscLin, pd.g.areaNuE0Tint);
    
  otherwise
    error('Unknown flux type')
end % switch

% Boundary edge matrices
globQOS = assembleMatEdgePhiIntPhiIntNu(pd.g, pd.g.markE0TbdrOS, refEdgePhiIntPhiInt, pd.g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu(pd.g, pd.g.markE0TbdrRA, refEdgePhiIntPhiInt, pd.g.areaNuE0TbdrRA);
globOL = assembleMatEdgePhiIntPhiIntFuncDiscNu(pd.g, pd.g.markE0TbdrL, refEdgePhiIntPhiIntPhiLin, zbDiscLin, pd.g.areaNuE0TbdrL);
globORA = assembleMatEdgePhiIntPhiIntFuncDiscNu(pd.g, pd.g.markE0TbdrRA, refEdgePhiIntPhiIntPhiLin, zbDiscLin, pd.g.areaNuE0TbdrRA);

% Derived system matrices
pd.sysW = blkdiag(pd.globM, pd.globM, pd.globM);
pd.linearTerms = [ sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1},  globQ{2} + globQOS{2} + globQRA{2} - globH{2}; ...
                   pd.gConst * (globG{1} + globU{1} - globO{1} - globOL{1} - globORA{1}), sparse(K*N,K*N), -globD; ...
                   pd.gConst * (globG{2} + globU{2} - globO{2} - globOL{2} - globORA{2}), globD, sparse(K*N,K*N) ];

% Slope limiting matrices
if ~isempty(pd.slopeLimList)
  globMTaylor = assembleMatElemPhiTaylorPhiTaylor(pd.g, N);
  pd.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(pd.g, N);
  pd.globMCorr = spdiags(1 ./ diag(globMTaylor), 0, K*N, K*N) * globMTaylor;
end % if

%% Assembly of time-independent global matrices corresponding to non-linear contributions.
% Element matrices
pd.globF = assembleMatElemDphiPerQuad(pd.g, refElemDphiPerQuad);

% Edge matrices
[pd.globRdiag, pd.globRoffdiag] = assembleMatEdgePhiNuPerQuad(pd.g, pd.g.markE0Tint, refEdgePhiIntPerQuad);
pd.globV = assembleMatEdgePhiPerQuad(pd.g, refEdgePhiIntPerQuad);

% Bottom-friction terms
if pd.isBottomFrictionVarying
  bottomFrictionDisc = projectFuncCont2DataDisc(pd.g, pd.bottomFrictionCont, 2*pd.p, pd.refElemPhiPhi);
  if pd.isBottomFrictionNonlinear
    refElemPhiPhiPerQuad = integrateRefElemPhiPhiPerQuad(N, pd.basesOnQuad);
    pd.globE = assembleMatElemPhiFuncDiscPerQuad(pd.g, refElemPhiPhiPerQuad, bottomFrictionDisc);
  else
    pd.globE = assembleMatElemPhiPhiFuncDisc(pd.g, refElemPhiPhiPhi, bottomFrictionDisc);
  end % if
else
  if pd.isBottomFrictionNonlinear
    refElemPhiPerQuad = integrateRefElemPhiPerQuad(N, pd.basesOnQuad);
    pd.globE = pd.bottomFrictionCoef * assembleMatElemPhiPhi(pd.g, refElemPhiPerQuad);
  else
    pd.globE = pd.bottomFrictionCoef * pd.globM;
  end % if
end % if

% Boundary matrices
if pd.g.numEbdrL > 0 % Land boundaries
  pd.globRL = assembleMatEdgePhiIntNuPerQuad(pd.g, pd.g.markE0TbdrL, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrL);
  switch pd.typeFluxL
    case 'reflected'
      error('not implemented')
    case 'natural'
      error('not implemented')
    case 'riemann'
      pd.globVL = assembleMatEdgePhiIntPerQuad(pd.g, pd.g.markE0TbdrL, refEdgePhiIntPerQuad, pd.g.areaE0TbdrL);
    otherwise
      error('Unknown flux type for land boundaries')
  end % switch
end % if

if pd.g.numEbdrRA > 0 % Radiation boundaries
  globRRA = assembleMatEdgePhiIntNuPerQuad(pd.g, pd.g.markE0TbdrRA, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrRA);
  pd.globRdiag = cellfun(@plus, pd.globRdiag, globRRA, 'UniformOutput', false);
end % if

if pd.g.numEbdrRI > 0 % River boundaries
  pd.globLRI = { sparse(K*N,1), sparse(K*N,1), sparse(K*N,1) };
  if ~pd.isRamping % TODO rhsRIAlg
    globRRI = assembleMatEdgePhiIntNuPerQuad(pd.g, pd.g.markE0TbdrRI, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrRI);
    for n = 1 : 3
      hRiv = pd.xiRiv(:,n) - pd.zbPerQuad{n};
      uHRiv = pd.uRiv(:,n) .* hRiv;
      vHRiv = pd.vRiv(:,n) .* hRiv;
      uvHRiv = uHRiv .* pd.vRiv(:,n);
      gHHRiv = pd.gConst * pd.xiRiv(:,n) .* ( 0.5 * pd.xiRiv(:,n) - pd.zbPerQuad{n} );
      pd.globLRI{1} = pd.globLRI{1} + globRRI{n,1} * uHRiv + globRRI{n,2} * vHRiv;
      pd.globLRI{2} = pd.globLRI{2} + globRRI{n,1} * (pd.uRiv(:,n) .* uHRiv + gHHRiv) + globRRI{n,2} * uvHRiv;
      pd.globLRI{3} = pd.globLRI{3} + globRRI{n,1} * uvHRiv + globRRI{n,2} * (pd.vRiv(:,n) .* vHRiv + gHHRiv);
    end % for
  end % if
end % if

if pd.g.numEbdrOS > 0 % Open sea boundaries
  pd.globROS = assembleMatEdgePhiIntNuPerQuad(pd.g, pd.g.markE0TbdrOS, refEdgePhiIntPerQuad, pd.g.areaNuE0TbdrOS);
  if pd.isRiemOS
    pd.globVOS = assembleMatEdgePhiIntPerQuad(pd.g, pd.g.markE0TbdrOS, refEdgePhiIntPerQuad, pd.g.areaE0TbdrOS);
  end % if
end % if

%% Assembly of rhs terms.
% Assemble Newtonian tide potential matrix
if pd.isTidalDomain
  for n = 1 : size(pd.forcingTidal, 3)
    for i = 1 : 2
      for j = 1 : 2
        pd.forcingTidal{i,j,n} = assembleMatElemPhiPhiFuncDiscConst(pd.g, pd.refElemPhiPhi, pd.forcingTidal{i,j,n});
      end % for
    end % for
  end % for
end % if

%% Variable timestepping preparation.
% TODO: variable time step
end % function
