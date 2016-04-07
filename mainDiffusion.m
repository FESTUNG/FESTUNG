% The main script for the diffusion equation as presented in @ref FRAK2015.
%
%===============================================================================
%> @file mainDiffusion.m
%>
%> @brief The main script for the diffusion equation as presented in @ref FRAK2015 .
%===============================================================================
%>
%> @brief The main script for the diffusion equation as presented in @ref FRAK2015.
%> 
%> This script builds a numerical solver to approximate solutions @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ of the diffusion equation
%> @f{align*}{
%> \mathbf{z}                                         &\;=\; - \nabla c                                  &&\text{in}~J\times\Omega\,,\\
%> \partial_t c  + \nabla\cdot (d\,\mathbf{z})        &\;=\; f                                           &&\text{in}~J\times\Omega\,,\\
%> c                                                  &\;=\; c_\mathrm{D}                                &&\text{on}~J\times{\partial\Omega}_{\mathrm{D}}\,,\\
%> \vec{z}\cdot\vec{\nu}                              &\;=\; g_\mathrm{N}                                &&\text{on}~J\times{\partial\Omega}_\mathrm{N}\,,\\
%> c                                                  &\;=\; c^0                                         &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The vector-valued quantity @f$\mathbf{z}@f$ was introduced as auxiliary unknown.  
%> The coefficients @f$d:J\times\Omega\rightarrow\mathbb{R}^+@f$ and @f$f:J\times\Omega\rightarrow \mathbb{R}@f$
%> may vary in time and space.  A detailed description is found in @ref FRAK2015.
%> 
%> This script can be used as a template for further modifications.
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function mainDiffusion()
more off % disable paging of output
tic % Start time measurement
%% Parameters.
hmax            = 1/8;        % maximum edge length of triangle
p               = 2;          % local polynomial degree
tEnd            = pi;         % end time
numSteps        = 20;         % number of time steps
isVisGrid       = true;       % visualization of grid
isVisSol        = true;       % visualization of solution
eta             = 1;          % penalty parameter (eta>0)
outputBasename  = 'solution'; % Basename of output files
outputTypes     = cellstr(['vtk';'tec']);
%% Parameter check.
diary([outputBasename '.log'])
assert(p >= 0 && p <= 4, 'Polynomial order must be zero to four.')
assert(hmax > 0        , 'Maximum edge length must be positive.' )
assert(numSteps > 0    , 'Number of time steps must be positive.')
%% Coefficients and boundary data.
c0     = @(x1,x2) sin(x1).*cos(x2);
dCont  = @(t,x1,x2) (x1<3/4&x1>1/4&x2<3/4&x2>1/4) + 0.01;
fCont  = @(t,x1,x2) 0.1*t*(x1==x1);
cDCont = @(t,x1,x2) sin(2*pi*x2 + t);
gNCont = @(t,x1,x2) x2;
%% Triangulation.
g = domainSquare(hmax);
if isVisGrid,  visualizeGrid(g);  end
%% Globally constant parameters.
K           = g.numT;                      % number of triangles
N           = nchoosek(p + 2, p);          % number of local DOFs
tau         = tEnd/numSteps;               % time step size
markE0Tint  = g.idE0T == 0;                % [K x 3] mark local edges that are interior
markE0TbdrN = g.idE0T == 1 | g.idE0T == 3; % [K x 3] mark local edges on the Neumann boundary
markE0TbdrD = ~(markE0Tint | markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary
g           = computeDerivedGridDataDiffusion(g, markE0Tint, markE0TbdrD, markE0TbdrN);
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K)
%% Lookup table for basis function.
computeBasesOnQuad(N);
%% Computation of matrices on the reference triangle.
hatM        = integrateRefElemPhiPhi(N);
hatG        = integrateRefElemDphiPhiPhi(N);
hatH        = integrateRefElemDphiPhi(N);
hatRdiag    = integrateRefEdgePhiIntPhiIntPhiInt(N);
hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(N);
hatSdiag    = integrateRefEdgePhiIntPhiInt(N);
hatSoffdiag = integrateRefEdgePhiIntPhiExt(N);
%% Assembly of time-independent global matrices.
globM  = assembleMatElemPhiPhi(g, hatM);
globH  = assembleMatElemDphiPhi(g, hatH);
globQ  = assembleMatEdgePhiPhiNu(g, markE0Tint, hatSdiag, hatSoffdiag, g.areaNuE0Tint);
globQN = assembleMatEdgePhiIntPhiIntNu(g, markE0TbdrN, hatSdiag);
globS  = eta * assembleMatEdgePhiPhi(g, markE0Tint, hatSdiag, hatSoffdiag);
globSD = eta * assembleMatEdgePhiIntPhiInt(g, markE0TbdrD, hatSdiag);
sysW = [ sparse(2*K*N,3*K*N) ; sparse(K*N,2*K*N), globM ];
%% Initial data.
cDisc = projectFuncCont2DataDisc(g, c0, 2*p, hatM);
sysY = [ zeros(2*K*N,1) ; reshape(cDisc', K*N, 1) ];
%% Time stepping.
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', tEnd, tau, numSteps)
for nStep = 1 : numSteps
  t = nStep * tau;
  %% L2-projections of algebraic coefficients.
  dDisc = projectFuncCont2DataDisc(g, @(x1,x2) dCont(t,x1,x2), 2*p, hatM);
  fDisc = projectFuncCont2DataDisc(g, @(x1,x2) fCont(t,x1,x2), 2*p, hatM);
  %% Assembly of time-dependent global matrices.
  globG = assembleMatElemDphiPhiFuncDisc(g, hatG, dDisc);
  globR = assembleMatEdgePhiPhiFuncDiscNu(g, markE0Tint, hatRdiag, hatRoffdiag, dDisc, g.areaNuE0Tint);
  %% Assembly of Dirichlet boundary contributions.
  globRD = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(g, markE0TbdrD, hatRdiag, dDisc, g.areaNuE0TbdrD);
  globJD = assembleVecEdgePhiIntFuncContNu(g, markE0TbdrD, @(x1,x2) cDCont(t,x1,x2), N, g.areaNuE0TbdrD);
  globKD = eta * assembleVecEdgePhiIntFuncCont(g, markE0TbdrD, @(x1,x2) cDCont(t,x1,x2), N);
  %% Assembly of Neumann boundary contributions.
  globKN = assembleVecEdgePhiIntFuncDiscIntFuncCont(g, markE0TbdrN, dDisc, @(x1,x2) gNCont(t,x1,x2), g.areaE0TbdrN);
  %% Assembly of the source contribution.
  globL = globM*reshape(fDisc', K*N, 1);
  %% Building and solving the system.
  sysA = [                        globM,              sparse(K*N,K*N), -globH{1}+globQ{1}+globQN{1};
                        sparse(K*N,K*N),                        globM, -globH{2}+globQ{2}+globQN{2};
           -globG{1}+globR{1}+globRD{1}, -globG{2}+globR{2}+globRD{2},                globS+globSD];
  sysV = [                   -globJD{1};                   -globJD{2};        globKD-globKN+globL];
  sysY = (sysW + tau*sysA) \ (sysW*sysY + tau*sysV);
  %% Visualization
  if isVisSol
    cDisc = reshape(sysY(2*K*N+1 : 3*K*N), N, K)';
    cLagr = projectDataDisc2DataLagr(cDisc);
    visualizeDataLagr(g, cLagr, 'c_h', outputBasename, nStep, outputTypes);
  end % if
end % for
fprintf('Total computation time: %g seconds.\n', toc);
diary off
end % function
