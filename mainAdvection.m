% The main script for the advection equation as presented in @ref RAWFK2016.
%
%===============================================================================
%> @file mainAdvection.m
%>
%> @brief The main script for the advection equation as presented in @ref RAWFK2016.
%===============================================================================
%>
%> @brief The main script for the adcvetion equation as presented in @ref RAWFK2016.
%> 
%> This script builds a numerical solver to approximate solutions @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ of the advection equation
%> @f{align*}{
%> \partial_t c  + \nabla\cdot (\mathbf{u}\,c)        &\;=\; f                                           &&\text{in}~J\times\Omega\,,\\
%> c                                                  &\;=\; c_\mathrm{D}                                &&\text{on}~J\times{\partial\Omega}_{\mathrm{in}}\,,\\
%> c                                                  &\;=\; c^0                                         &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The velocity @f$\mathbf{u}:J\times\Omega\rightarrow\mathbb{R}^2@f$ and 
%> right hand side@f$f:J\times\Omega\rightarrow \mathbb{R}@f$
%> may vary in time and space.  A detailed description is found in @ref RAWFK2016.
%> 
%> This script can be used as a template for further modifications.
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> Modified by Hennes Hajduk, 2016-04-06
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
function mainAdvection()
more off % disable paging of output
tic % Start time measurement
%% Parameters.
hmax        = 2^-6;              % maximum edge length of triangle
p           = 2;                 % local polynomial degree
ordRK       = min(p+1,3);        % order of Runge Kutta time stepper.
tEnd        = 2*pi;              % end time
numSteps    = 3142;              % number of time steps
isVisGrid   = false;             % visualization of grid
isVisSol    = true;              % visualization of solution
isSlopeLim  = true;              % slope limiting
typeSlopeLim = 'hierarch_vert';  % Type of slope limiter (linear, hierarch_vert, strict)
outputFrequency = 100;           % no visualization of every timestep
outputBasename  = ['solution_' typeSlopeLim]; % Basename of output files
outputTypes     = cellstr(['vtk';'tec']);
%% Parameter check.
diary([outputBasename '.log'])
assert(p >= 0 && p <= 4        , 'Polynomial order must be zero to four.'     )
assert(ordRK >= 1 && ordRK <= 3, 'Order of Runge Kutta must be zero to three.')
assert(hmax > 0                , 'Maximum edge length must be positive.'      )
assert(numSteps > 0            , 'Number of time steps must be positive.'     )
assert(~isSlopeLim || p > 0    , 'Slope limiting only available for p > 0.'   )
%% Triangulation.
g = domainSquare(hmax); % Alternative: g = domainPolygon([0 1 1 0], [0 0 1 1], hmax);
if isVisGrid,  visualizeGrid(g);  end
%% Globally constant parameters.
K           = g.numT;                      % number of triangles
N           = nchoosek(p + 2, p);          % number of local DOFs
tau         = tEnd/numSteps;               % time step size
markE0Tint  = g.idE0T == 0;                % [K x 3] mark local edges that are interior
markE0TbdrN = zeros(K, 3);                 % [K x 3] mark local edges on the Neumann boundary
markE0TbdrD = ~(markE0Tint | markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary
markV0TbdrD = ismember(g.V0T, g.V0E(g.E0T(markE0TbdrD),:)); % [K x 3] mark local vertices on the Dirichlet boundary
g           = computeDerivedGridDataAdvection(g, markE0TbdrD, markE0TbdrN);
%% Coefficients and boundary data (LeVeque's solid body rotation).
G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
c0Cont = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                    (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                   0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);
fCont = @(t,x1,x2) zeros(size(x1));
u1Cont = @(t,x1,x2) 0.5 - x2;
u2Cont = @(t,x1,x2) x1 - 0.5;
cDCont = @(t,x1,x2) zeros(size(x1));
gNCont = @(t,x1,x2) zeros(size(x1));
%% Lookup table for basis function.
computeBasesOnQuad(N);
if isSlopeLim
  computeTaylorBasesV0T(g, N);
end % if
%% Computation of matrices on the reference triangle.
hatM              = integrateRefElemPhiPhi(N);
hatG              = integrateRefElemDphiPhiPhi(N);
hatRdiagOnQuad    = integrateRefEdgePhiIntPhiIntPerQuad(N);
hatRoffdiagOnQuad = integrateRefEdgePhiIntPhiExtPerQuad(N);
%% Assembly of time-independent global matrices.
globM  = assembleMatElemPhiPhi(g, hatM);
if isSlopeLim
  globMTaylor     = assembleMatElemPhiTaylorPhiTaylor(g, N);
  globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(g, N);
  globMCorr       = spdiags(1./diag(globMTaylor), 0, K*N, K*N) * globMTaylor;
end % if
%% Initial data.
cDisc = projectFuncCont2DataDisc(g, c0Cont, 2*p+1, hatM);
if isSlopeLim
  cDV0T = computeFuncContV0T(g, @(x1, x2) cDCont(0, x1, x2));
  cDisc = applySlopeLimiterDisc(g, cDisc, markV0TbdrD, cDV0T, globM, globMDiscTaylor, typeSlopeLim);
end % if
fprintf('L2 error w.r.t. the initial condition: %g\n', computeL2Error(g, cDisc, c0Cont, 2*p));
%% visualization of inital condition.
if isVisSol
  cLagrange = projectDataDisc2DataLagr(cDisc);
  visualizeDataLagr(g, cLagrange, 'u_h', outputBasename, 0, outputTypes)
end
%% Time stepping.
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', tEnd, tau, numSteps)
for nStep = 1 : numSteps
  [t, omega] = rungeKuttaSSP(ordRK, tau, (nStep - 1) * tau);
  cDiscRK = cell(length(omega)+1, 1); cDiscRK{1} = reshape(cDisc', [K*N 1]);
  %% Perform Runge-Kutta steps
  for rkStep = 1 : length(omega)
    % L2 projections of Contebraic coefficients
    fDisc  = projectFuncCont2DataDisc(g, @(x1,x2) fCont(t(rkStep),x1,x2),  2*p, hatM);
    u1Disc = projectFuncCont2DataDisc(g, @(x1,x2) u1Cont(t(rkStep),x1,x2), 2*p, hatM);
    u2Disc = projectFuncCont2DataDisc(g, @(x1,x2) u2Cont(t(rkStep),x1,x2), 2*p, hatM);
    % Evaluate normal velocity in quadrature points of edges
    vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(g, @(x1,x2) u1Cont(t(rkStep),x1,x2), @(x1,x2) u2Cont(t(rkStep),x1,x2), 2*p+1); % veloc \dot \nu on quadratur points on edges
    % Assembly of time-dependent global matrices
    globG = assembleMatElemDphiPhiFuncDiscVec(g, hatG, u1Disc, u2Disc);
    globR = assembleMatEdgePhiPhiValUpwind(g, hatRdiagOnQuad, hatRoffdiagOnQuad, vNormalOnQuadEdge);
    % Assembly of Dirichlet boundary contributions
    globKD = assembleVecEdgePhiIntFuncContVal(g, markE0TbdrD, @(x1,x2) cDCont(t(rkStep),x1,x2), vNormalOnQuadEdge, N, g.areaE0TbdrD);
    % Assembly of Neumann boundary conditions
    globKN = assembleVecEdgePhiIntFuncCont(g, g.areaE0TbdrN, @(x1,x2) (gNCont(t(rkStep),x1,x2) <= 0) .* gNCont(t(rkStep),x1,x2), N);
    % Assembly of the source contribution
    globL = globM * reshape(fDisc', K*N, 1);
    % Building the system
    sysA = -globG{1} - globG{2} + globR;
    sysV = globL - globKD - globKN;
    % Computing the discrete time derivative
    cDiscDot = globM \ (sysV - sysA * cDiscRK{rkStep});
    % Apply slope limiting to time derivative
    if isSlopeLim
      cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot, [N K])', globM, globMDiscTaylor);
      cDiscDotTaylorLim = applySlopeLimiterTaylor(g, cDiscDotTaylor, markV0TbdrD, NaN(K,3), typeSlopeLim);
      cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
      cDiscDot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', globM, globMDiscTaylor)', [K*N 1]);
    end
    % Compute next step
    cDiscRK{rkStep + 1} = omega(rkStep) * cDiscRK{1} + (1 - omega(rkStep)) * (cDiscRK{rkStep} + tau * cDiscDot);
    % Limiting the solution
    if isSlopeLim
      cDV0T = computeFuncContV0T(g, @(x1, x2) cDCont(t(rkStep), x1, x2));
      cDiscRK{rkStep + 1} = reshape(applySlopeLimiterDisc(g, reshape(cDiscRK{rkStep + 1}, [N K])', markV0TbdrD, cDV0T, globM, globMDiscTaylor, typeSlopeLim)', [K*N 1]);
    end % if
  end % for
  cDisc = reshape(cDiscRK{end}, N, K)';
  %% visualization
  if isVisSol && mod(nStep, outputFrequency) == 0
    cLagrange = projectDataDisc2DataLagr(cDisc);
    visualizeDataLagr(g, cLagrange, 'u_h', outputBasename, nStep, outputTypes);
  end
end % for
if isVisSol
  cLagrange = projectDataDisc2DataLagr(cDisc);
  visualizeDataLagr(g, cLagrange, 'u_h', outputBasename, nStep, outputTypes);
end
fprintf('L2 error w.r.t. the initial condition: %g\n', computeL2Error(g, cDisc, c0Cont, 2*p));
fprintf('Total computation time: %g seconds.\n', toc);
diary off
end % function
