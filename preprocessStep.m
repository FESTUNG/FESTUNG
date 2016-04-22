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
function problemData = preprocessStep(problemData, nStep)
% Extract often used variables
K = problemData.K;
p = problemData.p;
N = problemData.N;
dt = problemData.dt;
t = problemData.t0 + nStep * dt;

% Determine time level at which continuous functions are to be evaluated
switch problemData.scheme
  case 'explicit'
    tRhs = t - dt;
  case 'semi_implicit'
    tRhs = t;
  otherwise
    error('Invalid time-stepping scheme.')  
end % switch

%% Source term contribution
if problemData.isRhsAvail
  % Project right hand side functions
  f0Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.f0Cont(x1,x2,tRhs), ...
            2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
  f1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.f1Cont(x1,x2,tRhs), ...
            2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
  f2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.f2Cont(x1,x2,tRhs), ...
            2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
  
  problemData.globL = cell(3,1);
  problemData.globL{1} = problemData.globM * reshape(f0Disc', K*N, 1);
  problemData.globL{2} = problemData.globM * reshape(f1Disc', K*N, 1);
  problemData.globL{3} = problemData.globM * reshape(f2Disc', K*N, 1);
elseif problemData.isTidalDomain
  error('not implemented')
else
  problemData.globL = { sparse(K*N,1); sparse(K*N,1); sparse(K*N,1) };
end % if

%% River boundary contributions
if problemData.isRiverBdr
  error('not implemented')
end % if

%% Compute water height on Open Sea boundaries.
qOrd1D = max(2*p, 1);  [Q, ~] = quadRule1D(qOrd1D);
heightOSPerQuad = cell(3,1);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  heightOSPerQuad{n} = problemData.xiOSCont(problemData.g.mapRef2Phy(1,Q1,Q2), problemData.g.mapRef2Phy(2,Q1,Q2), tRhs) - ...
                        problemData.zbPerQuad{n};
end % for

%% Create lookup tables for solution on quadrature points.
cDiscQ0T = cell(3,1); % cDisc in quadrature points of triangles
cDiscQ0E0Tint = cell(3,3); % cDisc in interior quad points of edges
cDiscQ0E0Text = cell(3,3,3); % cDisc in exterior quad points of edges
cDiscQ0E0TE0T = cell(3,3,3); % cDisc in quad points of edge on neighboring element

qOrd1D = 2*p + 1; [~, W] = quadRule1D(qOrd1D); numQuad1D = length(W);
qOrd2D = max(2*p,1); [~, ~, W] = quadRule2D(qOrd2D); numQuad2D = length(W);

for i = 1 : 3
  cDiscQ0T{i} = problemData.cDisc(:,:,i) * problemData.basesOnQuad.phi2D{qOrd2D}.';
  for nn = 1 : 3
    cDiscQ0E0Tint{i,nn} = problemData.cDisc(:,:,i) * problemData.basesOnQuad.phi1D{qOrd1D}(:,:,nn).';
    for np = 1 : 3
      cDiscQ0E0Text{i,nn,np} = problemData.cDisc(:,:,i) * problemData.basesOnQuad.thetaPhi1D{qOrd1D}(:,:,nn,np).';
      cDiscQ0E0TE0T{i,nn,np} = problemData.g.markE0TE0T{nn,np} * cDiscQ0E0Text{i,nn,np};
    end % for
  end % for
end % for

%%  Create lookup tables for eigenvalues for Lax-Friedrichs flux.
switch problemData.fluxType
  case 'Lax-Friedrichs'
    [lambdaE0TE0T, lambdaOSRiemE0T] = computeLaxFriedrichsCoefficientsSWE(problemData.g, problemData.gConst, cDiscQ0E0Tint, cDiscQ0E0TE0T, heightOSPerQuad, numQuad1D, problemData.averaging, problemData.isOSRiem);
  otherwise
    error('Unknown fluxType.')
end

%% Linearize lookup tables for solution and eigenvalues on quadrature points.
for i = 1 : 3
  cDiscQ0T{i} = reshape(cDiscQ0T{i}.', K * numQuad2D, 1);
  for nn = 1 : 3
    cDiscQ0E0Tint{i,nn} = reshape(cDiscQ0E0Tint{i,nn}.', K * numQuad1D, 1);
    for np = 1 : 3
      cDiscQ0E0Text{i,nn,np} = reshape(cDiscQ0E0Text{i,nn,np}.', K * numQuad1D, 1);
      cDiscQ0E0TE0T{i,nn,np} = reshape(cDiscQ0E0TE0T{i,nn,np}.', K * numQuad1D, 1);
    end % for
    lambdaE0TE0T{i,nn} = reshape(lambdaE0TE0T{i,nn}.', K * numQuad1D, 1);
  end % for
end % for

%% Evaluate non-linear terms on quadrature points.
% Non-linear terms in quadrature points of triangles
uuH = cDiscQ0T{2} .* cDiscQ0T{2} ./ cDiscQ0T{1};
uvH = cDiscQ0T{2} .* cDiscQ0T{3} ./ cDiscQ0T{1};
vvH = cDiscQ0T{3} .* cDiscQ0T{3} ./ cDiscQ0T{1};
gHH = 0.5 * problemData.gConst * (cDiscQ0T{1} .* cDiscQ0T{1});

problemData.nonlinearTerms = [ -problemData.globF{1} * (uuH + gHH) - problemData.globF{2} * uvH ; ...
                               -problemData.globF{1} * uvH - problemData.globF{2} * (vvH + gHH) ];
problemData.riemannTerms = sparse(3*K*N, 1);
                 
for nn = 1 : 3
  % Non-linear terms in exterior quadrature points of edges
  for np = 1 : 3
    uuH = cDiscQ0E0Text{2,nn,np} .* cDiscQ0E0Text{2,nn,np} ./ cDiscQ0E0Text{1,nn,np};
    uvH = cDiscQ0E0Text{2,nn,np} .* cDiscQ0E0Text{3,nn,np} ./ cDiscQ0E0Text{1,nn,np};
    vvH = cDiscQ0E0Text{3,nn,np} .* cDiscQ0E0Text{3,nn,np} ./ cDiscQ0E0Text{1,nn,np};
    gHH = 0.5 * problemData.gConst * (cDiscQ0E0Text{1,nn,np} .* cDiscQ0E0Text{1,nn,np});
    
    problemData.nonlinearTerms = problemData.nonlinearTerms + ...
                       [ problemData.globRoffdiag{nn,np,1} * (uuH + gHH) + problemData.globRoffdiag{nn,np,2} * uvH ; ...
                         problemData.globRoffdiag{nn,np,1} * uvH + problemData.globRoffdiag{nn,np,2} * (vvH + gHH) ];
                       
    problemData.riemannTerms = problemData.riemannTerms + ...
                       [ problemData.globV{nn,np} * (lambdaE0TE0T{nn,np} .* (cDiscQ0E0Tint{1,nn} - cDiscQ0E0TE0T{1,nn,np})) ; ...
                         problemData.globV{nn,np} * (lambdaE0TE0T{nn,np} .* (cDiscQ0E0Tint{2,nn} - cDiscQ0E0TE0T{2,nn,np})) ; ...
                         problemData.globV{nn,np} * (lambdaE0TE0T{nn,np} .* (cDiscQ0E0Tint{3,nn} - cDiscQ0E0TE0T{3,nn,np})) ];
  end % for
  
  % Non-linear terms in interior quadrature points of edges
  uuH = cDiscQ0E0Tint{2,nn} .* cDiscQ0E0Tint{2,nn} ./ cDiscQ0E0Tint{1,nn};
  uvH = cDiscQ0E0Tint{2,nn} .* cDiscQ0E0Tint{3,nn} ./ cDiscQ0E0Tint{1,nn};
  vvH = cDiscQ0E0Tint{3,nn} .* cDiscQ0E0Tint{3,nn} ./ cDiscQ0E0Tint{1,nn};
  gHH = 0.5 * problemData.gConst * (cDiscQ0E0Tint{1,nn} .* cDiscQ0E0Tint{1,nn});
  
  problemData.nonlinearTerms = problemData.nonlinearTerms + ...
                     [ problemData.globRdiag{nn,1} * (uuH + gHH) + problemData.globRdiag{nn,2} * uvH ; ...
                       problemData.globRdiag{nn,1} * uvH + problemData.globRdiag{nn,2} * (vvH + gHH) ];

  % Non-linear contributions of land boundaries
  problemData.nonlinearTerms = problemData.nonlinearTerms + ...
                     [ problemData.globRL{nn,1}; problemData.globRL{nn,2} ] * gHH;
                     
  % Non-linear contributions of Riemann solver on open sea edges
  if problemData.isOSRiem
     problemData.nonlinearTerms = problemData.nonlinearTerms + 0.5 * ...
                        [ problemData.globROS{nn,1} * (uuH + gHH) + problemData.globROS{nn,2} * uvH ; ...
                          problemData.globROS{nn,1} * uvH + problemData.globROS{nn,2} * (vvH + gHH) ];
  end % if
end % for

%% Evaluate bottom friction contributions.
if problemData.isBottomFrictionNonlinear
  normUoverH = sqrt(cDiscQ0T{2} .* cDiscQ0T{2} + cDiscQ0T{3} .* cDiscQ0T{3}) ./ ...
                              (cDiscQ0T{1} .* cDiscQ0T{1});
  problemData.bottomFrictionTerms = [ problemData.globE * (normUoverH .* cDiscQ0T{2}) ; ...
                                      problemData.globE * (normUoverH .* cDiscQ0T{3}) ];
else
  problemData.bottomFrictionTerms = [ problemData.globE * reshape(problemData.cDisc(:,:,2).', K*N, 1) ; ...
                                      problemData.globE * reshape(problemData.cDisc(:,:,3).', K*N, 1) ];
end % if

%% Other terms.
for i = 1 : 3
  for n = 1 : 3
    cDiscQ0E0Tint{i,n} = reshape(cDiscQ0E0Tint{i,n}, numQuad1D, K).';
  end % for
end % for
problemData.globUOSold = assembleGlobUOS(problemData.g, problemData.g.markE0TbdrOS, problemData.refEdgePhiIntPhiIntPerQuad, heightOSPerQuad, cDiscQ0E0Tint, problemData.g.areaE0TbdrOS);
problemData.globROSold = assembleGlobROS(problemData.g, problemData.g.markE0TbdrOS, heightOSPerQuad, N, problemData.gConst, problemData.basesOnQuad, problemData.g.areaNuE0TbdrOS);
if problemData.isOSRiem
  problemData.globVOSRiem = assembleGlobVOSRiem(problemData.g, problemData.g.markE0TbdrOS, heightOSPerQuad, problemData.refEdgePhiIntPhiIntPerQuad, lambdaOSRiemE0T, problemData.basesOnQuad, problemData.g.areaE0TbdrOS);
end % if
end % function
