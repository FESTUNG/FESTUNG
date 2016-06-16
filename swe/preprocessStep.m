% First step of the four-part algorithm in the main loop.

%===============================================================================
%> @file template/preprocessStep.m
%>
%> @brief First step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief First step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the number of
%> iterations provided by configureProblem in the parameter
%> <code>numSteps</code> is reached. These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed first in each loop iteration and is intended to
%> execute preprocessing operations, e.g., evaluate boundary conditions or
%> right hand side values, assemble time-dependent matrix blocks, etc.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
function pd = preprocessStep(pd, nStep)
% Extract often used variables
K = pd.K;
p = pd.p;
N = pd.N;

if pd.isVariableTimestep
  pd.dt = selectTimeStepSWE(pd.avgDiff(:,1), pd.avgDiff(:,2), pd.avgDepth, pd.gConst, pd.dt, nStep);
end % if

dt = pd.dt;
t = pd.t0 + nStep * dt;

% Determine time level at which continuous functions are to be evaluated
switch pd.schemeType
  case 'explicit'
    % TODO: Runge-Kutta
    tRhs = t - dt;
  case 'semi-implicit'
    tRhs = t;
  otherwise
    error('Invalid time-stepping scheme.')  
end % switch

%% Determine quadrature rules
qOrd1D = 2*p+1; [~, W] = quadRule1D(qOrd1D); numQuad1D = length(W);
qOrd2D = max(2*p,1); [~, ~, W] = quadRule2D(qOrd2D); numQuad2D = length(W);

%% Create lookup tables for solution on quadrature points.
cQ0T = cell(3,1); % cDisc in quadrature points of triangles
cQ0E0Tint = cell(3,3); % cDisc in interior quad points of edges
cQ0E0Text = cell(3,3,3); % cDisc in exterior quad points of edges
cQ0E0TE0T = cell(3,3,3); % cDisc in quad points of edge on neighboring element

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

hQ0T = cQ0T{1} - pd.zbQ0T; % water height (xi - zb) in quadrature points of triangles
% water height (xi - zb) in interior quad points of edges
hQ0E0Tint = cellfun(@minus, cQ0E0Tint(1,:).', pd.zbQ0E0Tint, 'UniformOutput', false);
% water height (xi - zb) in exterior quad points of edges
hQ0E0Text = cellfun(@minus, squeeze(cQ0E0Text(1,:,:)), pd.zbQ0E0Text, 'UniformOutput', false);
% water height (xi - zb) in exterior quad points of edge on neighboring element
hQ0E0TE0T = cellfun(@minus, squeeze(cQ0E0TE0T(1,:,:)), pd.zbQ0E0TE0T, 'UniformOutput', false);

markQ0E0TbdrL = cell(3,1);
markQ0E0TbdrOS = cell(3,1);
for n = 1 : 3
  markQ0E0TbdrL{n} = logical(kron(pd.g.markE0TbdrL(:,n), ones(numQuad1D,1)));
  markQ0E0TbdrOS{n} = logical(kron(pd.g.markE0TbdrOS(:,n), ones(numQuad1D,1)));
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
pd.tidalTerms = { sparse(K*N,K*N); sparse(K*N,K*N) };
if pd.isTidalDomain
  numFrequency = size(pd.forcingTidal, 3);
  for m = 1 : 2
    for n = 1 : numFrequency
      pd.tidalTerms{m} = pd.tidalTerms{m} + pd.forcingFrequency{1,n}(tRhs) * pd.forcingTidal{m,1,n} + ...
                                            pd.forcingFrequency{2,n}(tRhs) * pd.forcingTidal{m,2,n};
    end % for
    if pd.isRamp
      pd.tidalTerms{m} = pd.ramp(tRhs) * pd.tidalTerms{m};
    end % if
  end % for
end % if

%% Compute water height on Open Sea boundaries.
xiOSQ0E0Tint = cell(3,1);
if pd.g.numEbdrOS > 0
  if isfield(pd, 'xiOSCont')
    % Analytical function for open sea elevation given
    [Q, ~] = quadRule1D(max(2*p,1));
    for n = 1 : 3
      [Q1, Q2] = gammaMap(n, Q);
      xiOSQ0E0Tint{n} = pd.xiOSCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs);
      xiOSQ0E0Tint{n} = reshape(xiOSQ0E0Tint{n}.', K*numQuad1D,1);
    end % for
  elseif isfield(pd, 'xiFreqOS') && isfield(pd, 'xiAmpOS')
    % Open sea elevation data given
    % Since the open sea boundary condition is only used for non-linear
    % contributions we discretize it explicitly. Otherwise we would have
    % to make a distinction.
    numFrequency = size(pd.xiFreqOS, 2);
    xiOS = zeros(K, 1);
    for n = 1 : numFrequency
      xiOS = xiOS + pd.xiFreqOS{1,n}(tRhs) * pd.xiAmpOS{1,n} + pd.xiFreqOS{2,n}(tRhs) * pd.xiAmpOS{2,n};
    end % for
    xiOS = pd.ramp(tRhs) * kron(xiOS, ones(numQuad1D, 1));
    for n = 1 : 3
      xiOSQ0E0Tint{n} = xiOS;
    end % for
  else
    error('No open sea elevation given!')
  end % if
end % if

%% Compute river boundary values.
if pd.g.numEbdrRI > 0 && (pd.isRamp || pd.isRivCont)
  pd.globLRI = { sparse(K*N,1); sparse(K*N,1); sparse(K*N,1) };
  
  xiRivQ0E0T = cell(3,1);
  uRivQ0E0T = cell(3,1);
  vRivQ0E0T = cell(3,1);

  if pd.isRivCont
    [Q, ~] = quadRule1D(max(2*p,1));
    for n = 1 : 3
      [Q1, Q2] = gammaMap(n, Q);
      xiRivQ0E0T{n} = reshape(pd.xiRivCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
      uRivQ0E0T{n} = reshape(pd.uRivCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
      vRivQ0E0T{n} = reshape(pd.vRivCont(pd.g.mapRef2Phy(1,Q1,Q2), pd.g.mapRef2Phy(2,Q1,Q2), tRhs).', K*numQuad1D,1);
    end % for
  else
    for n = 1 : 3
      xiRivQ0E0T{n} = pd.ramp(tRhs) * pd.xiRivQ0E0T(:,n);
      uRivQ0E0T{n} = pd.ramp(tRhs) * pd.uRivQ0E0T(:,n);
      vRivQ0E0T{n} = pd.ramp(tRhs) * pd.vRivQ0E0T(:,n);
    end % for
  end % if
end % if

%% Evaluate bottom friction contributions.
if pd.isBottomFrictionNonlinear
  normUoverH = sqrt(cQ0T{2} .* cQ0T{2} + cQ0T{3} .* cQ0T{3}) ./ (hQ0T .* hQ0T);
  pd.bottomFrictionTerms = [ pd.globE * (normUoverH .* cQ0T{2}) ; pd.globE * (normUoverH .* cQ0T{3}) ];
else
  pd.bottomFrictionTerms = [ pd.globE * reshape(pd.cDisc(:,:,2).', K*N, 1) ; pd.globE * reshape(pd.cDisc(:,:,3).', K*N, 1) ];
end % if

%% Non-linear terms in quadrature points of triangles.
uuH = cQ0T{2} .* cQ0T{2} ./ hQ0T;
uvH = cQ0T{2} .* cQ0T{3} ./ hQ0T;
vvH = cQ0T{3} .* cQ0T{3} ./ hQ0T;
gHH = 0.5 * pd.gConst * (cQ0T{1} .* cQ0T{1});

pd.riemannTerms = sparse(3*K*N, 1);
pd.nonlinearTerms = [ -pd.globF{1} * (uuH + gHH) - pd.globF{2} * uvH ; ...
                      -pd.globF{1} * uvH - pd.globF{2} * (vvH + gHH) ];
                             
%% Non-linear terms in quadrature points of edges.
for nn = 1 : 3
  % Non-linear terms in exterior quadrature points of edges
  for np = 1 : 3
    uuH = cQ0E0Text{2,nn,np} .* cQ0E0Text{2,nn,np} ./ hQ0E0Text{nn,np};
    uvH = cQ0E0Text{2,nn,np} .* cQ0E0Text{3,nn,np} ./ hQ0E0Text{nn,np};
    vvH = cQ0E0Text{3,nn,np} .* cQ0E0Text{3,nn,np} ./ hQ0E0Text{nn,np};
    gHH = 0.5 * pd.gConst * (cQ0E0Text{1,nn,np} .* cQ0E0Text{1,nn,np});
    
    cAvgQ0E0T = computeAveragedVariablesQ0E0Tint(cQ0E0Tint(:,nn), cQ0E0TE0T(:,nn,np), hQ0E0Tint{nn}, hQ0E0TE0T{nn,np}, pd.averagingType);
    
    switch pd.typeFlux
      case 'Lax-Friedrichs'
        pd.nonlinearTerms = pd.nonlinearTerms + ...
          [ pd.globRoffdiag{nn,np,1} * (uuH + gHH) + pd.globRoffdiag{nn,np,2} * uvH ; ...
            pd.globRoffdiag{nn,np,1} * uvH + pd.globRoffdiag{nn,np,2} * (vvH + gHH) ];
          
        lambda = computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
        pd.riemannTerms = pd.riemannTerms + ...
          [ pd.globV{nn,np} * (lambda .* (cQ0E0Tint{1,nn} - cQ0E0TE0T{1,nn,np})) ; ...
            pd.globV{nn,np} * (lambda .* (cQ0E0Tint{2,nn} - cQ0E0TE0T{2,nn,np})) ; ...
            pd.globV{nn,np} * (lambda .* (cQ0E0Tint{3,nn} - cQ0E0TE0T{3,nn,np})) ];

      case 'Roe'
        error('not implemented')
        
      otherwise
        error('Invalid flux type for interior edges.')
    end % switch
  end % for
  
  % Non-linear terms in interior quadrature points of edges
  uuH = cQ0E0Tint{2,nn} .* cQ0E0Tint{2,nn} ./ hQ0E0Tint{nn};
  uvH = cQ0E0Tint{2,nn} .* cQ0E0Tint{3,nn} ./ hQ0E0Tint{nn};
  vvH = cQ0E0Tint{3,nn} .* cQ0E0Tint{3,nn} ./ hQ0E0Tint{nn};
  gHH = 0.5 * pd.gConst * (cQ0E0Tint{1,nn} .* cQ0E0Tint{1,nn});
  
  pd.nonlinearTerms = pd.nonlinearTerms + ...
    [ pd.globRdiag{nn,1} * (uuH + gHH) + pd.globRdiag{nn,2} * uvH ; ...
      pd.globRdiag{nn,1} * uvH + pd.globRdiag{nn,2} * (vvH + gHH) ];

  % Land boundary contributions
  if pd.g.numEbdrL > 0
    switch pd.typeBdrL
      case 'natural'
        pd.nonlinearTerms = pd.nonlinearTerms + [ pd.globRL{nn,1}; pd.globRL{nn,2} ] * gHH;

      case 'reflected'
        uHL = pd.g.nuE0Tsqr{nn,2} .* cQ0E0Tint{2,nn} - pd.g.nuE0Tprod{nn} .* cQ0E0Tint{3,nn};
        vHL = pd.g.nuE0Tsqr{nn,1} .* cQ0E0Tint{3,nn} - pd.g.nuE0Tprod{nn} .* cQ0E0Tint{2,nn};
        uuHL = uHL .* uHL ./ hQ0E0Tint{nn};
        uvHL = uHL .* vHL ./ hQ0E0Tint{nn};
        vvHL = vHL .* vHL ./ hQ0E0Tint{nn};
        % TODO: landFlux
        pd.nonlinearTerms = pd.nonlinearTerms + ...
          [ pd.globRdiag{nn,1} * (uuHL + gHH) + pd.globRdiag{nn,2} * uvHL ; ...
            pd.globRdiag{nn,1} * uvHL + pd.globRdiag{nn,2} * (vvHL + gHH) ];

      case 'riemann'
        uHriem = pd.g.nuE0TsqrDiff{nn} .* cQ0E0Tint{2,nn} - 2 * pd.g.nuE0Tprod{nn} .* cQ0E0Tint{3,nn};
        vHriem = -pd.g.nuE0TsqrDiff{nn} .* cQ0E0Tint{3,nn} - 2 * pd.g.nuE0Tprod{nn} .* cQ0E0Tint{2,nn};
        
        cQ0E0Triem = { [], uHriem, vHriem };
        cAvgQ0E0T = computeAveragedVariablesQ0E0Tland(cQ0E0Tint(:,nn), cQ0E0Triem, hQ0E0Tint{nn}, [], markQ0E0TbdrL{nn}, pd.averagingType);
        
        switch pd.typeFlux
          case 'Lax-Friedrichs'
            uuHriem = uHriem .* uHriem ./ hQ0E0Tint{nn};
            uvHriem = uHriem .* vHriem ./ hQ0E0Tint{nn};
            vvHriem = vHriem .* vHriem ./ hQ0E0Tint{nn};
            
            lambda = computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
            
            pd.nonlinearTerms = pd.nonlinearTerms + 0.5 * ...
              [ pd.globRL{nn,1} * (uuH + uuHriem + 2 * gHH) + pd.globRL{nn,2} * (uvH + uvHriem) + pd.globVL{nn} * (lambda .* (cQ0E0Tint{2,nn} - uHriem)); ...
                pd.globRL{nn,1} * (uvH + uvHriem) + pd.globRL{nn,2} * (vvH + vvHriem + 2 * gHH) + pd.globVL{nn} * (lambda .* (cQ0E0Tint{3,nn} - vHriem)) ];
                      % TODO: land flux
                      
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
  if pd.g.numEbdrRI > 0 && (pd.isRamp || pd.isRivCont)
    hRiv = xiRivQ0E0T{nn} - pd.zbQ0E0Tint{nn};
    uHRiv = uRivQ0E0T{nn} .* hRiv;
    vHRiv = vRivQ0E0T{nn} .* hRiv;
    uvHRiv = uRivQ0E0T{nn} .* vHRiv;
    gHHRiv = pd.gConst * xiRivQ0E0T{nn} .* (0.5 * xiRivQ0E0T{nn} - pd.zbQ0E0Tint{nn});
    
    pd.globLRI{1} = pd.globLRI{1} + pd.globRRI{nn,1} * uHRiv + pd.globRRI{nn,2} * vHRiv;
    pd.globLRI{2} = pd.globLRI{2} + pd.globRRI{nn,1} * (uRivQ0E0T{nn} .* uHRiv + gHHRiv) + pd.globRRI{nn,2} * uvHRiv;
    pd.globLRI{3} = pd.globLRI{3} + pd.globRRI{nn,1} * uvHRiv + pd.globRRI{nn,2} * (vRivQ0E0T{nn} .* vHRiv + gHHRiv);
  end % if
  
  % Open Sea boundary contributions
  if pd.g.numEbdrOS > 0
    hOSQ0E0T = xiOSQ0E0Tint{nn} - pd.zbQ0E0Tint{nn};
    validateattributes(hOSQ0E0T, {'numeric'}, {'>', 0})
    uuHOS = cQ0E0Tint{2,nn} .* cQ0E0Tint{2,nn} ./ hOSQ0E0T;
    uvHOS = cQ0E0Tint{2,nn} .* cQ0E0Tint{3,nn} ./ hOSQ0E0T;
    vvHOS = cQ0E0Tint{3,nn} .* cQ0E0Tint{3,nn} ./ hOSQ0E0T;
    gHHOS = pd.gConst * xiOSQ0E0Tint{nn} .* (0.5 * xiOSQ0E0Tint{nn} - pd.zbQ0E0Tint{nn});
    
    if pd.isRiemOS
      gHHOS = gHH + gHHOS - pd.gConst * cQ0E0Tint{1,nn} .* pd.zbQ0E0Tint{nn};
      
      pd.nonlinearTerms = pd.nonlinearTerms + 0.5 * ...
        [ pd.globROS{nn,1} * (uuH + uuHOS + gHHOS) + pd.globROS{nn,2} * (uvH + uvHOS) ; ...
          pd.globROS{nn,1} * (uvH + uvHOS) + pd.globROS{nn,2} * (vvH + vvHOS + gHHOS) ];
      
      cAvgQ0E0T = computeAveragedVariablesQ0E0Tos(cQ0E0Tint(:,nn), {}, hQ0E0Tint{nn}, hOSQ0E0T, markQ0E0TbdrOS{nn}, pd.averagingType);
        
      switch pd.typeFlux
        case 'Lax-Friedrichs'
          lambda = computeLaxFriedrichsCoefficient(cAvgQ0E0T, pd.g.nuQ0E0T(nn,:), pd.gConst);
          pd.riemannTerms(1:K*N) = pd.riemannTerms(1:K*N) + pd.globVOS{nn} * (lambda .* (hQ0E0Tint{nn} - hOSQ0E0T));
          
        case 'Roe'
          error('not implemented')
          
        otherwise
          error('Invalid flux type for land boundaries.')
      end % switch
    else
      pd.nonlinearTerms = pd.nonlinearTerms + ...
        [ pd.globROS{nn,1} * (uuHOS + gHHOS) + pd.globROS{nn,2} * uvHOS ; ...
          pd.globROS{nn,1} * uvHOS + pd.globROS{nn,2} * (vvHOS + gHHOS) ];
    end % if
  end % if
end % for
end % function
