% Sets all problem parameters and fills the problemData-struct with all basic 
% configuration options.

%===============================================================================
%> @file
%>
%> @brief Sets all problem parameters and fills the problemData-struct with
%>        all basic configuration options.
%===============================================================================
%>
%> @brief Sets all problem parameters and fills the problemData-struct with
%>        all basic configuration options.
%>
%> This routine is called before any other function for the problem.
%> It defines all problem parameters and should be the only file users have
%> to adapt to their needs.
%> 
%> Let @f$J = (0,t_\mathrm{end})@f$ be a finite time interval and @f$\Omega 
%> \subset \mathbb{R}^2@f$ a polygonally bounded domain with boundary 
%> @f$\partial\Omega@f$, subdivided into Dirichlet boundary parts 
%> @f$\partial\Omega_\mathrm{D}@f$ and Neumann boundary parts
%> @f$\partial\Omega_\mathrm{N}@f$.
%>
%> In this problem we consider the time-dependent Darcy equation
%> @f[
%> S_0 \partial_t h - \nabla \cdot (\mathsf{K} \nabla h) = \tilde{f}
%> @f]
%> describing water transport through fully saturated porous media, where @f$h@f$ 
%> is generally understood as the hydraulic head. 
%> The coefficient @f$S_0@f$ denotes the specific storativity of the porous medium
%> and @f$K@f$ the permeability. 
%> Division by @f$S_0@f$ and setting
%> @f[
%> \mathsf{D} := \mathsf{K} / S_0\,, \qquad
%> f := \tilde{f} / S_0
%> @f]
%> yields below formulation of the problem:
%>
%> We seek approximate solutions 
%> @f$c:\overline{J}\times\overline{\Omega}\rightarrow\mathbb{R}@f$ 
%> of the time-dependent diffusion equation
%> @f{align*}{
%> \mathbf{q}                                  &\;=\; - \nabla h   &&\text{in}~J\times\Omega\,,\\
%> \partial_t h  + \nabla\cdot (d\,\mathbf{q}) &\;=\; f            &&\text{in}~J\times\Omega\,,\\
%> c                                           &\;=\; h_\mathrm{D} &&\text{on}~J\times{\partial\Omega}_{\mathrm{D}}\,,\\
%> \vec{q}\cdot\vec{\nu}                       &\;=\; g_\mathrm{N} &&\text{on}~J\times{\partial\Omega}_\mathrm{N}\,,\\
%> c                                           &\;=\; h_0          &&\text{on}~\{0\}\times\Omega\,.
%> @f}
%> The vector-valued quantity @f$\mathbf{q}@f$ was introduced as auxiliary unknown.  
%> The coefficients @f$d:J\times\Omega\rightarrow\mathbb{R}^{2\times2}@f$ and 
%> @f$f:J\times\Omega\rightarrow \mathbb{R}@f$ may vary in time and space. 
%>
%> A detailed description can be found in @ref RRAFK2018.
%>
%> Please read the inline-comments in the code for the meaning of each
%> configuration option.
%>
%> @param  problemData  A (probably) empty struct with problem parameters.
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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
function problemData = configureProblem(problemData)
%% Parameters.
% Name of testcase
problemData = setdefault(problemData, 'testcase', 'showcase');

problemData.eta = 1; % penalty parameter (eta>0)

% Number of elements in x- and y-direction
problemData = setdefault(problemData, 'numElem', [32, 8]);

% Local polynomial approximation order (0 to 5)
problemData = setdefault(problemData, 'p', 1);

% Order of quadrature rule
problemData = setdefault(problemData, 'qOrd', 2*problemData.p + 1);
problemData = setdefault(problemData, 'qOrdMax', problemData.qOrd);

% Time stepping parameters
problemData = setdefault(problemData, 't0', 0);  % start time
problemData = setdefault(problemData, 'tEnd', 600);  % end time
problemData = setdefault(problemData, 'numSteps', ceil(problemData.tEnd/1));  % number of time steps

% Discard time derivative and compute stationary solution
problemData = setdefault(problemData, 'isStationary', false);  
problemData = setdefault(problemData, 'isCoupling', false);   
problemData = setdefault(problemData, 'isJumpCoupling', true);

% Visualization settings
problemData = setdefault(problemData, 'isVisGrid', false);  % visualization of grid
problemData = setdefault(problemData, 'isVisSol', true);  % visualization of solution
problemData = setdefault(problemData, 'outputFrequency', 100); % no visualization of every timestep
problemData = setdefault(problemData, 'outputBasename', ...  % Basename of output files
                         ['output' filesep problemData.problemName '_' problemData.testcase ]); 
problemData = setdefault(problemData, 'outputTypes', { 'vtk' });  % Type of visualization files ('vtk, 'tec')

%% Parameter check.
assert(problemData.p >= 0 && problemData.p <= 5, 'Polynomial order must be zero to five.')
assert(problemData.numSteps > 0, 'Number of time steps must be positive.')

%% Coefficients and boundary data.
isHotstart = false;

switch problemData.testcase    
  case 'showcase'
    isHotstart = false;
    hotstartFile = ['darcy_2dv' filesep 'showcase_p0_50x10.mat'];
    
    % width and height of computational domain
    zPMCont = @(x) -4 * ones(size(x));
    zBotCont = @(x) (cos((x-50)/30 * pi) + 1) .* (20 <= x & x <= 80);

    domainWidth = linspace(0, 100, problemData.numElem(1)+1);
    domainHeight = [zPMCont(domainWidth); zBotCont(domainWidth)];
    idDirichlet = 3; idNeumann = [1, 2, 4];
    
    k = 0.001;
    
    problemData.h0Cont = @(x,z) 5 * ones(size(x));
    problemData.q10Cont = @(x,z) zeros(size(x));
    problemData.q20Cont = @(x,z) zeros(size(x));
    
    problemData.hDCont = @(t,x,z) 5 * ones(size(x));
    problemData.gNCont = @(t,x,z) zeros(size(x));
    problemData.fCont = @(t,x,z) zeros(size(x));
    problemData.DCont = @(t,x,z) k * ones(size(x));

  case {'coupled_constXi', 'coupled_stationary', 'coupled_transient'}
    % Parameters
    aConst = 0;
    bConst = 0.005;
    betaConst = 0.3;
    etaConst = 0.003;
    rhoConst = 0.08;
    tauConst = 0.08;
    kConst = 0.01;
    kappaConst = 1;
    lambdaConst = 0.07;
    muConst = 0.07;
    
    % Width and height of computational domain
    zPMCont = @(x) -5 * ones(size(x));
    zBotCont = @(x) aConst + bConst * x;
    
    domainWidth = linspace(0, 100, problemData.numElem(1)+1);
    domainHeight = [zPMCont(domainWidth); zBotCont(domainWidth)];
    idDirichlet = [1, 2, 3, 4]; idNeumann = -1;
    
    % Analytical solution
    switch problemData.testcase
      case 'coupled_constXi'
        xiCont = @(t,x) 5 * ones(size(x));

        dxXiCont = @(t,x) zeros(size(x));
        dxdxXiCont = @(t,x) zeros(size(x));
        dtXiCont = @(t,x) zeros(size(x));
      case 'coupled_stationary'
        xiCont = @(t,x) 5 + etaConst * sin(rhoConst * x);

        dxXiCont = @(t,x) etaConst * rhoConst * cos(rhoConst * x);
        dxdxXiCont = @(t,x) -etaConst * rhoConst * rhoConst * sin(rhoConst * x);
        dtXiCont = @(t,x) zeros(size(x));
      case 'coupled_transient'
        xiCont = @(t,x) 5 + etaConst * sin(rhoConst * x + tauConst * t);

        dxXiCont = @(t,x) etaConst * rhoConst * cos(rhoConst * x + tauConst * t);
        dxdxXiCont = @(t,x) -etaConst * rhoConst * rhoConst * sin(rhoConst * x + tauConst * t);
        dtXiCont = @(t,x) etaConst * tauConst * cos(rhoConst * x + tauConst * t);
    end % switch
    
    switch problemData.testcase
      case {'coupled_constXi', 'coupled_stationary'}
        omegaCont = @(t,x) kappaConst * cos(lambdaConst * x);
        dtOmegaCont = @(t,x) zeros(size(x));
        dxOmegaCont = @(t,x) -lambdaConst * kappaConst * sin(lambdaConst * x);
        dxdxOmegaCont = @(t,x) - lambdaConst * lambdaConst * kappaConst * cos(lambdaConst * x);
      case {'coupled_transient'}
        omegaCont = @(t,x) kappaConst * cos(lambdaConst * x + muConst * t);
        dtOmegaCont = @(t,x) -muConst * kappaConst * sin(lambdaConst * x + muConst * t);
        dxOmegaCont = @(t,x) -lambdaConst * kappaConst * sin(lambdaConst * x + muConst * t);
        dxdxOmegaCont = @(t,x) -lambdaConst * lambdaConst * kappaConst * cos(lambdaConst * x + muConst * t);
    end % switch
    
    problemData.hCont = @(t,x,z) xiCont(t,x) + (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* omegaCont(t,x);
    problemData.q1Cont = @(t,x,z) -dxXiCont(t,x) ...
      + betaConst * bConst * cos(betaConst * zBotCont(x)) .* omegaCont(t,x) ...
      - (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* dxOmegaCont(t,x);
    problemData.q2Cont = @(t,x,z) -betaConst * cos(betaConst * z) .* omegaCont(t, x);
    
    % Diffusion matrix. Can be specified as:
    % - 2x2 cell
    % - 2-element cell, will be interpreted as diagonal matrix
    % - scalar function, will be interpreted as diagonal matrix
    problemData.DCont = @(t,x,z) kConst * ones(size(x));
    
    % Derivatives
    dThCont = @(t,x,z) dtXiCont(t,x) + (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* dtOmegaCont(t,x);
    dXhCont = @(t,x,z) -problemData.q1Cont(t,x,z);
    dZhCont = @(t,x,z) -problemData.q2Cont(t,x,z);
    dXdXhCont = @(t,x,z) -dxdxXiCont(t,x) ...
      + betaConst * betaConst * bConst * bConst * sin(betaConst * zBotCont(x)) .* omegaCont(t,x) ...
      - 2 * betaConst * bConst * cos(betaConst * zBotCont(x)) .* dxOmegaCont(t,x) ...
      + (sin(betaConst * z) - sin(betaConst * zBotCont(x))) .* dxdxOmegaCont(t,x);
    dZdZhCont = @(t,x,z) -betaConst * betaConst * sin(betaConst * z) .* omegaCont(t, x);
    dXZDCont = cellfun(@(c) @(t,x,z) c * ones(size(x)), {0, 0; 0, 0}, 'UniformOutput', false);
    
    % Boundary conditions
    problemData.hDCont = problemData.hCont;
    problemData.gNCont = @(t,x,z) zeros(size(x));
    
    % Right hand side
    problemData.fCont = @(t,x,z) dThCont(t,x,z) - ...
                          dXZDCont{1}(t,x,z) .* dXhCont(t,x,z)  - problemData.DCont(t,x,z) .* dXdXhCont(t,x,z) - ...
                          dXZDCont{2}(t,x,z) .* dZhCont(t,x,z)  - problemData.DCont(t,x,z) .* dZdZhCont(t,x,z);
                   
  case 'convergence'
    % width and height of computational domain
    domainWidth = [0, 10];
    domainHeight = [0, 10];
    idDirichlet = [1, 2, 3, 4]; idNeumann = -1;
    % Analytical solution
    problemData.hCont = @(t,x,z) cos(x + t) .* cos(z + t);
    problemData.q1Cont = @(t,x,z) sin(x + t) .* cos(z + t);
    problemData.q2Cont = @(t,x,z) cos(x + t) .* sin(z + t);
    % Diffusion matrix
    problemData.DCont = { @(t,x,z) exp(z/5) , @(t,x,z) 0.5 * ones(size(x)) ; ...
                          @(t,x,z) 1/3 * ones(size(x)), @(t,x,z) exp(x/5) };
    % Derivatives
    dThCont = @(t,x,z) -sin(x + t) .* cos(z + t) - cos(x + t) .* sin(z + t);
    dXhCont = @(t,x,z) -sin(x + t) .* cos(z + t);
    dZhCont = @(t,x,z) -cos(x + t) .* sin(z + t);
    dXdXhCont = @(t,x,z) -cos(x + t) .* cos(z + t);
    dZdZhCont = @(t,x,z) -cos(x + t) .* cos(z + t);
    dXdZhCont = @(t,x,z) sin(x + t) .* sin(z + t);
    dXZDCont = { @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)); ...
                 @(t,x,z) zeros(size(x)), @(t,x,z) zeros(size(x)) };
    % Boundary conditions
    problemData.hDCont = problemData.hCont;
    problemData.gNCont = @(t,x,z) zeros(size(x));
    % Right hand side
    problemData.fCont = @(t,x,z) dThCont(t,x,z) - ...
                          dXZDCont{1,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.DCont{1,1}(t,x,z) .* dXdXhCont(t,x,z) - ...
                          dXZDCont{1,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.DCont{1,2}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZDCont{2,1}(t,x,z) .* dXhCont(t,x,z)  - problemData.DCont{2,1}(t,x,z) .* dXdZhCont(t,x,z) - ...
                          dXZDCont{2,2}(t,x,z) .* dZhCont(t,x,z)  - problemData.DCont{2,2}(t,x,z) .* dZdZhCont(t,x,z);
                        
  case 'convergence_stationary'
    problemData.isStationary = true;
    % width and height of computational domain
    domainWidth = [0, 1];
    domainHeight = [0, 1];
    idDirichlet = [1, 2, 3, 4]; idNeumann = -1;
    % Analytical solution
    cx = 7 / domainWidth(2); cz = 7 / domainHeight(2);
    dx = 1 / domainWidth(2); dz = 1 / domainHeight(2);
    problemData.hCont = @(t,x,z) cos(cx * x) .* cos(cz * z);
    problemData.q1Cont = @(t,x,z) cx * sin(cx * x) .* cos(cz * z);
    problemData.q2Cont = @(t,x,z) cz * cos(cx * x) .* sin(cz * z);
    % Diffusion matrix
    problemData.DCont = @(t,x,z) exp(dx * x + dz * z);
    % Derivatives
    dXhCont = @(t,x,z) -cx * sin(cx * x) .* cos(cz * z);
    dZhCont = @(t,x,z) -cz * cos(cx * x) .* sin(cz * z);
    dXdXhCont = @(t,x,z) -cx * cx * cos(cx * x) .* cos(cz * z);
    dZdZhCont = @(t,x,z) -cz * cz * cos(cx * x) .* cos(cz * z);
    dXDCont = @(t,x,z) dx * exp(dx * x + dz * z);
    dZDCont = @(t,x,z) dz * exp(dx * x + dz * z);
    % Boundary conditions
    problemData.hDCont = problemData.hCont;
    problemData.gNCont = @(t,x,z) zeros(size(x));
    % Right hand side
    problemData.fCont = @(t,x,z) -dXDCont(t,x,z) .* dXhCont(t,x,z) - dZDCont(t,x,z) .* dZhCont(t,x,z) ...
      - problemData.DCont(t,x,z) .* (dXdXhCont(t,x,z) + dZdZhCont(t,x,z));

  otherwise
    error('Unknown testcase "%s".\n', problemData.testcase);
end % switch

problemData = setdefault(problemData, 'isHotstart', isHotstart);
if problemData.isHotstart, problemData = setdefault(problemData, 'hotstartFile', hotstartFile); end

problemData = setdefault(problemData, 'idNeumann', idNeumann);
problemData = setdefault(problemData, 'idDirichlet', idDirichlet);

%% Domain and triangulation.
problemData.generateGrid = @(numElem) domainRectTrap(domainWidth, domainHeight, numElem);

% Boundary parts (0 = int, 1 = bot, 2 = right, 3 = top, 4 = left)
checkMultipleIds = @(idE0T, ids) logical(sum(bsxfun(@eq, idE0T, reshape(ids, 1, 1, length(ids))), 3));

problemData.generateMarkE0Tint = @(g) g.idE0T == 0;
problemData.generateMarkE0TbdrCoupling = @(g) g.idE0T == 3;
problemData.generateMarkE0TbdrN = @(g) checkMultipleIds(g.idE0T, problemData.idNeumann) & ~(problemData.isCoupling & g.idE0T == 3);
problemData.generateMarkE0TbdrD = @(g) checkMultipleIds(g.idE0T, problemData.idDirichlet) & ~(problemData.isCoupling & g.idE0T == 3);
end % function