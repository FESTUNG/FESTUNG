% First step of the three-part algorithm in the iterateSubSteps loop of each
% step in the main loop.

%===============================================================================
%> @file
%>
%> @brief First step of the three-part algorithm in the iterateSubSteps loop of 
%>        each step in the main loop.
%===============================================================================
%>
%> @brief First step of the three-part algorithm in the iterateSubSteps loop of 
%>        each step in the main loop.
%>
%> The iterateSubSteps loop repeatedly executes three steps until the number of
%> substep iterations equals the order of the underlying Runge-Kutta method.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%> 
%> This routine is executed first in each substep loop iteration and is intended
%> to execute preprocessing operations, e.g., evaluate boundary conditions or
%> right hand side values, assemble time-dependent matrix blocks, etc.
%>
%> @param  pd           A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval pd           The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
function pd = preprocessSubStep(pd, ~, nSubStep)
% Extract often used variables
K = pd.K;
p = pd.p;
N = pd.N;

tRhs = pd.tLvls(nSubStep);

%% Determine quadrature rules
qOrd1D = 2*p+1; [~, W] = quadRule1D(qOrd1D); numQuad1D = length(W);
qOrd2D = max(2*p,1); [~, ~, W] = quadRule2D(qOrd2D); numQuad2D = length(W);

%% Create lookup tables for solution on quadrature points.
cQ0T = cell(5,1); % cDisc in quadrature points of triangles
cQ0E0Tint = cell(5,3); % cDisc in interior quad points of edges
cQ0E0Text = cell(5,3,3); % cDisc in exterior quad points of edges
cQ0E0TE0T = cell(5,3,3); % cDisc in quad points of edge on neighboring element

for i = 1 : 3
  cQ0T{i} = reshape(pd.basesOnQuad.phi2D{qOrd2D} * pd.cDisc(:,:,i).', numQuad2D * K, 1);
  for nn = 1 : 3
    cQ0E0Tint{i,nn} = reshape(pd.basesOnQuad.phi1D{qOrd1D}(:,:,nn) * pd.cDisc(:,:,i).', numQuad1D * K, 1);
    for np = 1 : 3
      cDiscThetaPhi = pd.basesOnQuad.thetaPhi1D{qOrd1D}(:,:,nn,np) * pd.cDisc(:,:,i).';
      cQ0E0Text{i,nn,np} = reshape(cDiscThetaPhi, numQuad1D * K, 1);
      cQ0E0TE0T{i,nn,np} = reshape(cDiscThetaPhi * pd.g.markE0TE0T{nn,np}.', numQuad1D * K, 1);
    end % for
  end % for
end % for

cQ0T{4} = reshape(pd.basesOnQuad.phi2D{qOrd2D} * pd.rampInput(pd.t) * pd.xiDisc(:,:).', numQuad2D * K, 1);
cQ0T{5} = reshape(pd.basesOnQuad.phi2D{qOrd2D} * pd.zbDisc(:,:).', numQuad2D * K, 1);

% cQ0T{1} = cQ0T{4} - cQ0T{5};

for nn = 1 : 3
  cQ0E0Tint{4,nn} = reshape(pd.basesOnQuad.phi1D{qOrd1D}(:,:,nn) * pd.rampInput(pd.t) * pd.xiDisc(:,:).', numQuad1D * K, 1);
  cQ0E0Tint{5,nn} = reshape(pd.basesOnQuad.phi1D{qOrd1D}(:,:,nn) * pd.zbDisc(:,:).', numQuad1D * K, 1);
  
%   cQ0E0Tint{1,nn} = cQ0E0Tint{4,nn} - cQ0E0Tint{5,nn};
  
  for np = 1 : 3
    cDiscThetaPhi = pd.basesOnQuad.thetaPhi1D{qOrd1D}(:,:,nn,np) * pd.rampInput(pd.t) * pd.xiDisc(:,:).';
    cQ0E0Text{4,nn,np} = reshape(cDiscThetaPhi, numQuad1D * K, 1);
    cQ0E0TE0T{4,nn,np} = reshape(cDiscThetaPhi * pd.g.markE0TE0T{nn,np}.', numQuad1D * K, 1);
		cDiscThetaPhi = pd.basesOnQuad.thetaPhi1D{qOrd1D}(:,:,nn,np) * pd.zbDisc(:,:).';
    cQ0E0Text{5,nn,np} = reshape(cDiscThetaPhi, numQuad1D * K, 1);
    cQ0E0TE0T{5,nn,np} = reshape(cDiscThetaPhi * pd.g.markE0TE0T{nn,np}.', numQuad1D * K, 1);

%     cQ0E0Text{1,nn,np} = cQ0E0Text{4,nn,np} - cQ0E0Text{5,nn,np};
%     cQ0E0TE0T{1,nn,np} = cQ0E0TE0T{4,nn,np} - cQ0E0TE0T{5,nn,np};
    
  end % for
end % for

%% ramping boundary condition for bathymtry
zbD = cell(3,1);
zbD{1} = -pd.constValue + pd.rampInput(pd.t) * (pd.zbQ0E0Tint{1} + pd.constValue);
zbD{2} = -pd.constValue + pd.rampInput(pd.t) * (pd.zbQ0E0Tint{2} + pd.constValue);
zbD{3} = -pd.constValue + pd.rampInput(pd.t) * (pd.zbQ0E0Tint{3} + pd.constValue);

markQ0E0TbdrL = cell(3,1);
markQ0E0TbdrRI = cell(3,1);
markQ0E0TbdrOS = cell(3,1);
markQ0E0TbdrF  = cell(3,1);
for n = 1 : 3
  markQ0E0TbdrL{n} = logical(kron(pd.g.markE0TbdrL(:,n), ones(numQuad1D,1)));
  markQ0E0TbdrRI{n} = logical(kron(pd.g.markE0TbdrRI(:,n), ones(numQuad1D,1)));
  markQ0E0TbdrOS{n} = logical(kron(pd.g.markE0TbdrOS(:,n), ones(numQuad1D,1)));
  markQ0E0TbdrF {n} = logical(kron(pd.g.markE0TbdrF (:,n), ones(numQuad1D,1)));
end % for

%% Right hand side contribution
pd.globL = { sparse(K*N,1); sparse(K*N,1); sparse(K*N,1) };
if pd.isRhsAvail
  f0Disc = projectFuncCont2DataDisc(pd.g, @(x1,x2) pd.f0Cont(x1,x2,tRhs), 2*p, pd.refElemPhiPhi, pd.basesOnQuad);
  f1Disc = projectFuncCont2DataDisc(pd.g, @(x1,x2) pd.f1Cont(x1,x2,tRhs), 2*p, pd.refElemPhiPhi, pd.basesOnQuad);
  f2Disc = projectFuncCont2DataDisc(pd.g, @(x1,x2) pd.f2Cont(x1,x2,tRhs), 2*p, pd.refElemPhiPhi, pd.basesOnQuad);
  
  pd.globL{1} = pd.globM * reshape(f0Disc.', K*N, 1);
  pd.globL{2} = pd.globM * reshape(f1Disc.', K*N, 1);
  pd.globL{3} = pd.globM * reshape(f2Disc.', K*N, 1);
end % if

%% Tidal potential contribution
pd.tidalTerms = { sparse(K*N,K*max(N,3)); sparse(K*N,K*max(N,3)) };
if pd.isTidalDomain
  numFrequency = size(pd.forcingTidal, 3);
  for m = 1 : 2
    for n = 1 : numFrequency
      pd.tidalTerms{m} = pd.tidalTerms{m} + pd.forcingFrequency{1,n}(tRhs) * pd.forcingTidal{m,1,n} + ...
                                            pd.forcingFrequency{2,n}(tRhs) * pd.forcingTidal{m,2,n};
    end % for
    if pd.isRamp
      pd.tidalTerms{m} = pd.ramp(tRhs/86400) * pd.tidalTerms{m};
    end % if
  end % for
end % if

%% Compute water height on Open Sea boundaries.
xiOSQ0E0T = cell(3,1);
if pd.g.numEbdrOS > 0
  if pd.isOSCont
    % Analytical function for open sea elevation given
    [Q, ~] = quadRule1D(max(2*p,1));
    for n = 1 : 3
      [Q1, Q2] = gammaMap(n, Q);
      xiOSQ0E0T{n} = pd.rampInput(pd.t) * pd.xiOSCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs);
      xiOSQ0E0T{n} = reshape(xiOSQ0E0T{n}.', K*numQuad1D,1);
    end % for
    if any(ismember(pd.slopeLimList, 'elevation'))
      pd.xiV0Tos = pd.g.markV0TbdrOS .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiOSCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)));
    end % if
  elseif isfield(pd, 'xiFreqOS') && isfield(pd, 'xiAmpOS')
    % Open sea elevation data given
    % Since the open sea boundary condition is only used for non-linear
    % contributions we discretize it explicitly. Otherwise we would have
    % to make a distinction.
    numFrequency = size(pd.xiFreqOS, 2);
    xiOS = zeros(K, 3);
    xiE0Tos = zeros(K, 3);
    for n = 1 : numFrequency
      xiOS = xiOS + pd.xiFreqOS{1,n}(tRhs) * pd.xiAmpOS{1,n} + pd.xiFreqOS{2,n}(tRhs) * pd.xiAmpOS{2,n};
      if any(ismember(pd.slopeLimList, 'elevation'))
        xiE0Tos = xiE0Tos + pd.xiFreqOS{1,n}(getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)) * pd.xiAmpOS{1,n} + ...
                            pd.xiFreqOS{2,n}(getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)) * pd.xiAmpOS{2,n};
      end % if
    end % for
    if any(ismember(pd.slopeLimList, 'elevation'))
      % To determine the vertex values for all triangles we first compute
      % the vertex values for triangles that have a boundary edge of open sea
      % type. We average the values from both sides of the vertex if the 
      % triangle has more than one boundary edge.
      xiV0T = [ sum(xiE0Tos(:,[2,3]),2) ./ sum(pd.g.markE0TbdrOS(:,[2,3]),2), ...
                sum(xiE0Tos(:,[1,3]),2) ./ sum(pd.g.markE0TbdrOS(:,[1,3]),2), ...
                sum(xiE0Tos(:,[1,2]),2) ./ sum(pd.g.markE0TbdrOS(:,[1,2]),2) ];

      % the vertex values for triangles with boundary edges and (possibly) 
      % NaN if it is not a boundary vertex
      xiV0T = xiV0T(pd.g.markV0TbdrOS);

      xiV = pd.vertInd2VertIndUniqueOS * setNaN2Zero(xiV0T);
      
      % xiV contains the sum of the vertex values of all triangles and is
      % divided by the number of contributing elements in dataVCountOS
      xiV = xiV ./ pd.xiVCountOS;
      pd.xiV0Tos = xiV(pd.g.V0T);
    end % if
    for n = 1 : 3
      xiOSQ0E0T{n} = pd.ramp(tRhs/86400) * pd.rampInput(pd.t) * kron(xiOS(:,n), ones(numQuad1D, 1));
    end % for
  else
    error('No open sea elevation given.')
  end % if
end % if

%% Compute river boundary values.
if pd.g.numEbdrRI > 0 && (pd.isRamp || pd.isRivCont || pd.isRiemRiv)
  pd.globLRI = { sparse(K*N,1); sparse(K*N,1); sparse(K*N,1) };
  
  xiRivQ0E0T = cell(3,1);
  uRivQ0E0T = cell(3,1);
  vRivQ0E0T = cell(3,1);

  if pd.isRivCont
    [Q, ~] = quadRule1D(max(2*p,1));
    for n = 1 : 3
      [Q1, Q2] = gammaMap(n, Q);
      xiRivQ0E0T{n} = pd.rampInput(pd.t) * reshape(pd.xiRivCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
      uRivQ0E0T{n} = reshape(pd.uRivCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
      vRivQ0E0T{n} = reshape(pd.vRivCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
    end % for
    if any(ismember(pd.slopeLimList, 'momentum'))
      xiV0T = computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)));
      if any(ismember(pd.slopeLimList, 'elevation'))
        pd.xiV0Triv = pd.g.markV0TbdrRI .* xiV0T;
      end % if
      hV0T = xiV0T - pd.zbV0T;
      pd.uHV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.uRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt))) .* hV0T;
      pd.vHV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.vRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt))) .* hV0T;
    elseif any(ismember(pd.slopeLimList, 'elevation'))
      pd.xiV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)));
    end % if
  else
    for n = 1 : 3
      xiRivQ0E0T{n} = pd.rampInput(pd.t) * pd.ramp(tRhs/86400) * pd.xiRivQ0E0T(:,n);
      uRivQ0E0T{n} = pd.ramp(tRhs/86400) * pd.uRivQ0E0T(:,n);
      vRivQ0E0T{n} = pd.ramp(tRhs/86400) * pd.vRivQ0E0T(:,n);
    end % for
  end % if
end % if

%% Compute flow boundary values.
pd.globLF = { sparse(K*N,1); sparse(K*N,1); sparse(K*N,1) };
if pd.g.numEbdrF > 0
  
  uHFQ0E0T = cell(3,1);
  vHFQ0E0T = cell(3,1);

  if pd.isFlowCont
    [Q, ~] = quadRule1D(max(2*p,1));
    for n = 1 : 3
      [Q1, Q2] = gammaMap(n, Q);
      uHFQ0E0T{n} = reshape(pd.uHFCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
      vHFQ0E0T{n} = reshape(pd.vHFCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
    end % for
%     if any(ismember(pd.slopeLimList, 'momentum'))
%       xiV0T = computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)));
%       if any(ismember(pd.slopeLimList, 'elevation'))
%         pd.xiV0Triv = pd.g.markV0TbdrRI .* xiV0T;
%       end % if
%       hV0T = xiV0T - pd.zbV0T;
%       pd.uHV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.uRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt))) .* hV0T;
%       pd.vHV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.vRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt))) .* hV0T;
%     elseif any(ismember(pd.slopeLimList, 'elevation'))
%       pd.xiV0Triv = pd.g.markV0TbdrRI .* computeFuncContV0T(pd.g, @(x1, x2) pd.xiRivCont(x1, x2, getdefault(pd.tLvls, nSubStep+1, pd.t + pd.dt)));
%     end % if
  else
    for n = 1 : 3
      uHFQ0E0T{n} = pd.ramp(tRhs/86400) * pd.uHFQ0E0T(:,n);
      vHFQ0E0T{n} = pd.ramp(tRhs/86400) * pd.vHFQ0E0T(:,n);
    end % for
  end % if
end % if

%% Evaluate nudging term
% zbDiff = pd.observationVertices .* (-pd.depth - pd.zbV);
% pd.nudgingTerm = 0.0 * -pd.rampInput(pd.t) * pd.matElemPhiPhiLagr * reshape(zbDiff(pd.g.V0T).', K*N, 1);

%% Evaluate bottom friction contributions.
if pd.isBottomFrictionNonlinear
  normUoverH = sqrt(cQ0T{2} .* cQ0T{2} + cQ0T{3} .* cQ0T{3}) ./ (cQ0T{1} .* cQ0T{1});
  pd.bottomFrictionTerms = [ pd.globE * (normUoverH .* cQ0T{2}) ; pd.globE * (normUoverH .* cQ0T{3}) ];
else
  pd.bottomFrictionTerms = [ pd.globE * reshape(pd.cDisc(:,:,2).', K*N, 1) ; pd.globE * reshape(pd.cDisc(:,:,3).', K*N, 1) ];
end % if

%% Non-linear terms in quadrature points of triangles.
uuH = cQ0T{2} .* cQ0T{2} ./ cQ0T{1};
uvH = cQ0T{2} .* cQ0T{3} ./ cQ0T{1};
vvH = cQ0T{3} .* cQ0T{3} ./ cQ0T{1};
gEE = 0.5 * pd.gConst * (cQ0T{4} .* cQ0T{4});

pd.riemannTerms = sparse(3*K*N, 1);
pd.nonlinearTerms = [ -pd.globF{1} * (uuH+gEE) - pd.globF{2} * uvH ; ...
                      -pd.globF{1} * uvH - pd.globF{2} * (vvH+gEE) ];

%% assemble gravitational term
globG = assembleMatElemPhiPhiDfuncDisc(pd.g, pd.refElemPhiPhiDphi, pd.zbDisc);
globU = assembleMatElemDphiPhiFuncDisc(pd.g, pd.refElemDphiPhiPhi, pd.zbDisc);
% globO can only be computed this way in case of Lax-Friedrichs flux
globO   = assembleMatEdgePhiPhiFuncDiscIntNu(pd.g, pd.g.markE0Tint, pd.refEdgePhiIntPhiIntPhi, pd.refEdgePhiIntPhiExtPhi, pd.zbDisc, pd.g.areaNuE0Tint);
globOL  = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(pd.g, pd.g.markE0TbdrL, pd.refEdgePhiIntPhiIntPhi, pd.zbDisc, pd.g.areaNuE0TbdrL);
globORA = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(pd.g, pd.g.markE0TbdrRA, pd.refEdgePhiIntPhiIntPhi, pd.zbDisc, pd.g.areaNuE0TbdrRA);
globOF  = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(pd.g, pd.g.markE0TbdrF, pd.refEdgePhiIntPhiIntPhi, pd.zbDisc, pd.g.areaNuE0TbdrF);

pd.gravityTerm = pd.gConst * [ globG{1} + globU{1} - globO{1} - globOL{1} - globORA{1} - globOF{1} ; ...
                               globG{2} + globU{2} - globO{2} - globOL{2} - globORA{2} - globOF{2} ] * pd.rampInput(pd.t) * reshape(pd.xiDisc.', K*N, 1);

if pd.isCoupling
  pd.massFluxQ0E0T = zeros(K, 3, numQuad1D);
  pd.uHDisc = pd.cDisc(:,:,2);
  pd.vHDisc = pd.cDisc(:,:,3);
end % if

%% Non-linear terms in quadrature points of edges.
for nn = 1 : 3
  if pd.isCoupling
    pd.massFluxQ0E0T(:,nn,:) = 0.5 * bsxfun(@times, permute( reshape( ( cQ0E0Tint{2,nn} + sum(cat(2,cQ0E0TE0T{2,nn,:}),2) ) .* pd.g.nuQ0E0T{nn,1} + ...
                                                                      ( cQ0E0Tint{3,nn} + sum(cat(2,cQ0E0TE0T{3,nn,:}),2) ) .* pd.g.nuQ0E0T{nn,2}, ...
                                                                      [numQuad1D, K, 1] ), [2 3 1] ), pd.g.markE0Tint(:,nn) ) ...
                               + bsxfun(@times, permute( reshape( cQ0E0Tint{2,nn} .* pd.g.nuQ0E0T{nn,1} + cQ0E0Tint{3,nn} .* pd.g.nuQ0E0T{nn,2}, ...
                                                                  [numQuad1D, K, 1] ), [2 3 1] ), pd.g.markE0TbdrRA(:,nn) );
  end % if
  % Non-linear terms in exterior quadrature points of edges
  for np = 1 : 3
    uuH = cQ0E0Text{2,nn,np} .* cQ0E0Text{2,nn,np} ./ cQ0E0Text{1,nn,np};
    uvH = cQ0E0Text{2,nn,np} .* cQ0E0Text{3,nn,np} ./ cQ0E0Text{1,nn,np};
    vvH = cQ0E0Text{3,nn,np} .* cQ0E0Text{3,nn,np} ./ cQ0E0Text{1,nn,np};
    gEE = 0.5 * pd.gConst * (cQ0E0Text{4,nn,np} .* cQ0E0Text{4,nn,np});
		
    cAvgQ0E0T = pd.computeAveragedVariablesQ0E0Tint(cQ0E0Tint(1:3,nn), cQ0E0TE0T(1:3,nn,np), cQ0E0Tint{1,nn}, cQ0E0TE0T{1,nn,np}, pd.averagingType);
    
    switch pd.typeFlux
      case 'Lax-Friedrichs'
        pd.nonlinearTerms = pd.nonlinearTerms + ...
          [ pd.globRoffdiag{nn,np,1} * (uuH+gEE) + pd.globRoffdiag{nn,np,2} * uvH ; ...
            pd.globRoffdiag{nn,np,1} * uvH + pd.globRoffdiag{nn,np,2} * (vvH+gEE) ];
          
        lambda = pd.computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
        pd.riemannTerms = pd.riemannTerms + ...
          [ pd.globV{nn,np} * (lambda .* (cQ0E0Tint{1,nn} - cQ0E0TE0T{1,nn,np})) ; ...
            pd.globV{nn,np} * (lambda .* (cQ0E0Tint{2,nn} - cQ0E0TE0T{2,nn,np})) ; ...
            pd.globV{nn,np} * (lambda .* (cQ0E0Tint{3,nn} - cQ0E0TE0T{3,nn,np})) ];

        if pd.isCoupling
          if isfield(pd.g, 'markE0T')
            pd.massFluxQ0E0T(:,nn,:) = pd.massFluxQ0E0T(:,nn,:) + 0.5 * bsxfun(@times, permute( reshape( lambda .* (cQ0E0Tint{1,nn} ...
                                                                          - cQ0E0TE0T{1,nn,np}), [numQuad1D, K, 1] ), [2 3 1] ), pd.g.markE0T{nn,np} );
          else
            pd.massFluxQ0E0T(:,nn,:) = pd.massFluxQ0E0T(:,nn,:) + 0.5 * bsxfun(@times, permute( reshape( lambda .* (cQ0E0Tint{1,nn} ...
                                                            - cQ0E0TE0T{1,nn,np}), [numQuad1D, K, 1] ), [2 3 1] ), pd.g.markE0TE0T{nn,np} * ones(K,1) );
          end % if
        end % if
      case 'Roe'
        error('not implemented')
        
      otherwise
        error('Invalid flux type for interior edges.')
    end % switch
  end % for
  
  % Non-linear terms in interior quadrature points of edges
  uuH = cQ0E0Tint{2,nn} .* cQ0E0Tint{2,nn} ./ cQ0E0Tint{1,nn};
  uvH = cQ0E0Tint{2,nn} .* cQ0E0Tint{3,nn} ./ cQ0E0Tint{1,nn};
  vvH = cQ0E0Tint{3,nn} .* cQ0E0Tint{3,nn} ./ cQ0E0Tint{1,nn};
  gEE = 0.5 * pd.gConst * (cQ0E0Tint{4,nn} .* cQ0E0Tint{4,nn});
  
  pd.nonlinearTerms = pd.nonlinearTerms + [ pd.globRdiag{nn,1} * (uuH+gEE) + pd.globRdiag{nn,2} * uvH ; ...
                                            pd.globRdiag{nn,1} * uvH + pd.globRdiag{nn,2} * (vvH+gEE) ];

  % Land boundary contributions
  if pd.g.numEbdrL > 0
    switch pd.typeBdrL
      case 'natural'
        pd.nonlinearTerms = pd.nonlinearTerms + [ pd.globRL{nn,1}; pd.globRL{nn,2} ] * gEE;
        
        % for coupled problems massFluxQ0E0T is zero on land boundary edges
        % with natural discretization

      case 'reflected'
        uHL = pd.g.nuE0Tsqr{nn,2} .* cQ0E0Tint{2,nn} - pd.g.nuE0Tprod{nn} .* cQ0E0Tint{3,nn};
        vHL = pd.g.nuE0Tsqr{nn,1} .* cQ0E0Tint{3,nn} - pd.g.nuE0Tprod{nn} .* cQ0E0Tint{2,nn};
        uuHL = uHL .* uHL ./ cQ0E0Tint{1,nn};
        uvHL = uHL .* vHL ./ cQ0E0Tint{1,nn};
        vvHL = vHL .* vHL ./ cQ0E0Tint{1,nn};
        
        pd.nonlinearTerms = pd.nonlinearTerms + ...
          [ pd.globRL{nn,1} * (uuHL+gEE) + pd.globRL{nn,2} * uvHL ; ...
            pd.globRL{nn,1} * uvHL + pd.globRL{nn,2} * (vvHL+gEE) ];
        
        % for coupled problems massFluxQ0E0T is zero on land boundary edges
        % with reflected discretization

      case 'riemann'
        uHriem = pd.g.nuE0TsqrDiff{nn} .* cQ0E0Tint{2,nn} - 2 * pd.g.nuE0Tprod{nn} .* cQ0E0Tint{3,nn};
        vHriem = -pd.g.nuE0TsqrDiff{nn} .* cQ0E0Tint{3,nn} - 2 * pd.g.nuE0Tprod{nn} .* cQ0E0Tint{2,nn};
        
        cQ0E0Triem = { [], uHriem, vHriem };
        cAvgQ0E0T = pd.computeAveragedVariablesQ0E0Tland(cQ0E0Tint(1:3,nn), cQ0E0Triem, cQ0E0Tint{1,nn}, [], markQ0E0TbdrL{nn}, pd.averagingType);
        
        switch pd.typeFlux
          case 'Lax-Friedrichs'
            uuHriem = uHriem .* uHriem ./ cQ0E0Tint{1,nn};
            uvHriem = uHriem .* vHriem ./ cQ0E0Tint{1,nn};
            vvHriem = vHriem .* vHriem ./ cQ0E0Tint{1,nn};
            
            lambda = pd.computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
            
            pd.nonlinearTerms = pd.nonlinearTerms + 0.5 * ...
              [ pd.globRL{nn,1} * (uuH + uuHriem + 2*gEE) + pd.globRL{nn,2} * (uvH + uvHriem) + pd.globVL{nn} * (lambda.*(cQ0E0Tint{2,nn}-uHriem)); ...
                pd.globRL{nn,1} * (uvH + uvHriem) + pd.globRL{nn,2} * (vvH + vvHriem + 2*gEE) + pd.globVL{nn} * (lambda.*(cQ0E0Tint{3,nn}-vHriem)) ];
            
            % for coupled problems massFluxQ0E0T is zero on land boundary edges
            % with Lax-Friedrichs-flux discretization
            
          case 'Roe'
            error('not implemented')
                                          
          otherwise
            error('Invalid flux type for land boundaries.')
        end % switch
        
      otherwise
        error('Invalid type for land boundary treatment.')
    end % switch
  end % if

  % River boundary contributions
  if pd.g.numEbdrRI > 0 && (pd.isRamp || pd.isRivCont || pd.isRiemRiv)
    
    hRiv = xiRivQ0E0T{nn} - zbD{nn};
    uHRiv = uRivQ0E0T{nn} .* hRiv;
    vHRiv = vRivQ0E0T{nn} .* hRiv;
    
    if pd.isRiemRiv
      cAvgQ0E0T = pd.computeAveragedVariablesQ0E0Triv(cQ0E0Tint(1:3,nn), { [], uHRiv, vHRiv }, cQ0E0Tint{1,nn}, hRiv, markQ0E0TbdrRI{nn}, pd.averagingType);
      switch pd.typeFlux
        case 'Lax-Friedrichs'
          lambda = pd.computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
          
          uuHRiv = uuH + uRivQ0E0T{nn} .* uHRiv;
          uvHRiv = uvH + uRivQ0E0T{nn} .* vHRiv;
          vvHRiv = vvH + vRivQ0E0T{nn} .* vHRiv;
          gEERiv = gEE + pd.gConst * (-cQ0E0Tint{4,nn} .* cQ0E0Tint{5,nn} + xiRivQ0E0T{nn} .* (0.5 * xiRivQ0E0T{nn} - zbD{nn}));
          
          uRiv = cQ0E0Tint{2,nn} + uHRiv;
          vRiv = cQ0E0Tint{3,nn} + vHRiv;
          jump = lambda .* (cQ0E0Tint{1,nn} - hRiv);

          pd.globLRI{1} = pd.globLRI{1} + 0.5 * ( pd.globRRI{nn,1} * uRiv + pd.globRRI{nn,2} * vRiv + pd.globVRI{nn} * jump );
          pd.globLRI{2} = pd.globLRI{2} + 0.5 * ( pd.globRRI{nn,1} * (uuHRiv+gEERiv) + pd.globRRI{nn,2} * uvHRiv ... 
                                                + pd.globVRI{nn} * (lambda .* (cQ0E0Tint{2,nn} - uHRiv)) );
          pd.globLRI{3} = pd.globLRI{3} + 0.5 * ( pd.globRRI{nn,1} * uvHRiv + pd.globRRI{nn,2} * (vvHRiv+gEERiv) ...
                                                + pd.globVRI{nn} * (lambda .* (cQ0E0Tint{3,nn} - vHRiv)) );
          
          if pd.isCoupling
            pd.massFluxQ0E0T(:,nn,:) = pd.massFluxQ0E0T(:,nn,:) + bsxfun(@times, 0.5 * permute( reshape( uRiv .* pd.g.nuQ0E0T{nn,1} + ...
                                              vRiv .* pd.g.nuQ0E0T{nn,2} + jump, [numQuad1D, K, 1] ), [2 3 1] ), pd.g.markE0TbdrRI(:,nn) );
          end % if
          
        case 'Roe'
          error('not implemented')
          
        otherwise
          error('Invalid flux type for river boundaries.')
      end % switch

    else
      uuHRiv = uRivQ0E0T{nn} .* uHRiv;
      uvHRiv = uRivQ0E0T{nn} .* vHRiv;
      vvHRiv = vRivQ0E0T{nn} .* vHRiv;
      gEERiv = pd.gConst * xiRivQ0E0T{nn} .* (0.5 * xiRivQ0E0T{nn} - zbD{nn});

      pd.globLRI{1} = pd.globLRI{1} + pd.globRRI{nn,1} * uHRiv + pd.globRRI{nn,2} * vHRiv;
      pd.globLRI{2} = pd.globLRI{2} + pd.globRRI{nn,1} * (uuHRiv+gEERiv) + pd.globRRI{nn,2} * uvHRiv;
      pd.globLRI{3} = pd.globLRI{3} + pd.globRRI{nn,1} * uvHRiv + pd.globRRI{nn,2} * (vvHRiv+gEERiv);
      
      if pd.isCoupling % squeeze needed since river boundary flux is sparse
        pd.massFluxQ0E0T(:,nn,:) = squeeze(pd.massFluxQ0E0T(:,nn,:)) + bsxfun(@times, reshape( uHRiv .* pd.g.nuQ0E0T{nn,1} ...
                                                   + vHRiv .* pd.g.nuQ0E0T{nn,2}, [numQuad1D, K] ).', pd.g.markE0TbdrRI(:,nn) );
      end % if
    end % if
  end % if
  
  % Open Sea boundary contributions
  if pd.g.numEbdrOS > 0
    hOSQ0E0T = xiOSQ0E0T{nn} - cQ0E0Tint{5,nn};
    validateattributes(hOSQ0E0T, {'numeric'}, {'>', 0})
    uuHOS = cQ0E0Tint{2,nn} .* cQ0E0Tint{2,nn} ./ hOSQ0E0T;
    uvHOS = cQ0E0Tint{2,nn} .* cQ0E0Tint{3,nn} ./ hOSQ0E0T;
    vvHOS = cQ0E0Tint{3,nn} .* cQ0E0Tint{3,nn} ./ hOSQ0E0T;
    gEEOS = pd.gConst * xiOSQ0E0T{nn} .* (0.5 * xiOSQ0E0T{nn} - cQ0E0Tint{5,nn});
    
    if pd.isRiemOS
      cAvgQ0E0T = pd.computeAveragedVariablesQ0E0Tos(cQ0E0Tint(1:3,nn), {}, cQ0E0Tint{1,nn}, hOSQ0E0T, markQ0E0TbdrOS{nn}, pd.averagingType);
      
      switch pd.typeFlux
        case 'Lax-Friedrichs'
          gEEOS = gEE + gEEOS - pd.gConst * cQ0E0Tint{4,nn} .* cQ0E0Tint{5,nn};
      
          pd.nonlinearTerms = pd.nonlinearTerms + 0.5 * [ pd.globROS{nn,1} * (uuH + uuHOS + gEEOS) + pd.globROS{nn,2} * (uvH + uvHOS); ...
                                                          pd.globROS{nn,1} * (uvH + uvHOS) + pd.globROS{nn,2} * (vvH + vvHOS + gEEOS) ];
          
          lambda = pd.computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
          pd.riemannTerms(1:K*N) = pd.riemannTerms(1:K*N) + 0.5 * pd.globVOS{nn} * (lambda .* (cQ0E0Tint{1,nn} - hOSQ0E0T));
          
          if pd.isCoupling
            pd.massFluxQ0E0T(:,nn,:) = pd.massFluxQ0E0T(:,nn,:) + bsxfun(@times, permute( reshape( cQ0E0Tint{2,nn} .* pd.g.nuQ0E0T{nn,1} + ...
                                          cQ0E0Tint{3,nn} .* pd.g.nuQ0E0T{nn,2} + 0.5 * lambda .* (cQ0E0Tint{1,nn} - hOSQ0E0T), [numQuad1D, K, 1] ), ...
                                          [2 3 1] ), pd.g.markE0TbdrOS(:,nn) );
          end % if
          
        case 'Roe'
          error('not implemented')
          
        otherwise
          error('Invalid flux type for open sea boundaries.')
      end % switch
    else
			gEEOS = pd.gConst * xiOSQ0E0T{nn} .* (0.5 * xiOSQ0E0T{nn} - cQ0E0Tint{5,nn});
      pd.nonlinearTerms = pd.nonlinearTerms + ...
        [ pd.globROS{nn,1} * (uuHOS+gEEOS) + pd.globROS{nn,2} * uvHOS ; ...
          pd.globROS{nn,1} * uvHOS + pd.globROS{nn,2} * (vvHOS+gEEOS) ];
      
      if pd.isCoupling
        pd.massFluxQ0E0T(:,nn,:) = pd.massFluxQ0E0T(:,nn,:) + bsxfun(@times, permute( reshape( cQ0E0Tint{2,nn} .* pd.g.nuQ0E0T{nn,1} + ...
                                          cQ0E0Tint{3,nn} .* pd.g.nuQ0E0T{nn,2}, [numQuad1D, K, 1] ), [2 3 1] ), pd.g.markE0TbdrOS(:,nn) );
      end % if
    end % if
  end % if
  
  % Flow boundary contributions
  if pd.g.numEbdrF > 0
    
    hF = cQ0E0Tint{4,nn} - zbD{nn};
    uHF = uHFQ0E0T{nn};
    vHF = vHFQ0E0T{nn};
    
    if pd.isRiemFlow
      cAvgQ0E0T = pd.computeAveragedVariablesQ0E0Tflow(cQ0E0Tint(1:3,nn), { [], uHF, vHF }, cQ0E0Tint{1,nn}, {}, markQ0E0TbdrF{nn}, pd.averagingType);
      switch pd.typeFlux
        case 'Lax-Friedrichs'
          lambda = pd.computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
          
          uuHF = uuH + uHF .* uHF ./ hF;
          uvHF = uvH + uHF .* vHF ./ hF;
          vvHF = vvH + vHF .* vHF ./ hF;
          
          pd.globLF{1} = pd.globLF{1} + 0.5 * ( pd.globRF{nn,1} * (uHF  + cQ0E0Tint{2,nn}) + pd.globRF{nn,2} * (vHF  + cQ0E0Tint{3,nn}) + pd.globVF{nn} * (lambda .* (cQ0E0Tint{1,nn} - hF)) );
          pd.globLF{2} = pd.globLF{2} + 0.5 * ( pd.globRF{nn,1} * (uuHF + 2 * gEE) + pd.globRF{nn,2} * uvHF + pd.globVF{nn} * (lambda .* (cQ0E0Tint{2,nn} - uHF)));
          pd.globLF{3} = pd.globLF{3} + 0.5 * ( pd.globRF{nn,1} * uvHF + pd.globRF{nn,2} * (vvHF + 2 * gEE) + pd.globVF{nn} * (lambda .* (cQ0E0Tint{3,nn} - vHF)));
        case 'Roe'
          error('not implemented')
          
        otherwise
          error('Invalid flux type for river boundaries.')
      end % switch
    else
      uuHF = uHF .* uHF ./ hF;
      uvHF = uHF .* vHF ./ hF;
      vvHF = vHF .* vHF ./ hF;

      pd.globLF{1} = pd.globLF{1} + pd.globRF{nn,1} * uHF + pd.globRF{nn,2} * vHF;
      pd.globLF{2} = pd.globLF{2} + pd.globRF{nn,1} * (uuHF + gEE) + pd.globRF{nn,2} * uvHF;
      pd.globLF{3} = pd.globLF{3} + pd.globRF{nn,1} * uvHF + pd.globRF{nn,2} * (vvHF + gEE);
    end % if
  end % if
  
end % for

if pd.isCoupling && (~pd.isRamp && ~pd.isRivCont && ~pd.isRiemRiv)
  pd.massFluxQ0E0T = pd.massFluxQ0E0T + pd.massFluxQ0E0TRiv;
end % if

end % function
