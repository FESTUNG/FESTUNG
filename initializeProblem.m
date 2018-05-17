% Fills the problem's data structures with initial data.

%===============================================================================
%> @file
%>
%> @brief Fills the problem's data structures with initial data.
%===============================================================================
%>
%> @brief Fills the problem's data structures with initial data.
%>
%> This routine is called after swe_2dv/preprocessProblem.m
%>
%> Before entering the main loop of the solution algorithm, this routine
%> fills the problem's data structures with initial data.
%>
%> It projects the initial condition and initializes the solution vectors.
%> If <tt>problemData.isHotstart</tt> is set to true, it loads the current
%> state from the file specified in <tt>problemData.hotstart</tt>.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by swe_2dv/configureProblem.m and 
%>                      swe_2dv/preprocessProblem.m. @f$[\text{struct}]@f$
%>
%> @retval problemData  The input struct enriched with initial data.
%>                      @f$[\text{struct}]@f$
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
function problemData = initializeProblem(problemData)
problemData.isFinished = false;

% Vector of unknowns (H, U, W) for each Runge-Kutta stage
problemData.cDiscRK = cell(length(rungeKuttaExplicit(problemData.ordRK, 0, 0)), 3);

%% Load state
if problemData.isHotstart
  hotstart = load(problemData.hotstartFile);
  assert(all(isfield(hotstart, {'t', 'hDisc', 'u1Disc', 'u2Disc'})), 'Hotstart file must contain t, hDisc, u1Disc, u2Disc')
  assert(problemData.g.g1D.numT == size(hotstart.hDisc, 1), 'Number of 1D elements in hotstart file does not match!')
  assert(all(problemData.g.numT == [size(hotstart.u1Disc, 1), size(hotstart.u2Disc, 1)]), 'Number of 2D elements in hotstart file does not match!')
  fprintf('Loaded hotstart data from "%s" at time level t=%g.\n', problemData.hotstartFile, hotstart.t);
end % if

%% Initial height.
if problemData.isHotstart
  barN = min(problemData.barN, size(hotstart.hDisc, 2));
  problemData.cDiscRK{end, 1} = zeros(problemData.g.g1D.numT, problemData.barN);
  problemData.cDiscRK{end, 1}(:, 1:barN) = hotstart.hDisc(:, 1:barN);
else
  problemData.cDiscRK{end, 1} = projectFuncCont2DataDisc1D(problemData.g.g1D, problemData.h0Cont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D);
end % if

%% Mesh adaptivity and assembly of time-independent global matrices.
problemData = problemData.fn_adaptFreeSurface(problemData, true);

%% Computation of bathymetry gradient.
dZbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 2);
dXbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 1) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 1);
dXzBot = problemData.g.g1D.markT2DT * ( problemData.gConst * (dZbot1D ./ dXbot1D) );
problemData.globLzBot = kron(dXzBot, eye(problemData.N, 1));

%% Initial velocities.
if problemData.isHotstart
  N = min(problemData.N, size(hotstart.u1Disc, 2));
  problemData.cDiscRK{end, 2} = zeros(problemData.g.numT, problemData.N);
  problemData.cDiscRK{end, 2}(:, 1:N) = hotstart.u1Disc(:, 1:N);
  problemData.cDiscRK{end, 3} = zeros(problemData.g.numT, problemData.N);
  problemData.cDiscRK{end, 3}(:, 1:N) = hotstart.u2Disc(:, 1:N);
else
  problemData.cDiscRK{end, 2} = projectFuncCont2DataDiscTetra(problemData.g, problemData.u10Cont, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
  % problemData.cDiscRK{end, 3} = projectFuncCont2DataDiscTetra(problemData.g, @(x,z) problemData.u2Cont(0,x,z), problemData.qOrd, problemData.globM, problemData.basesOnQuad2D);
  problemData.cDiscRK{end, 3} = zeros(size(problemData.cDiscRK{end, 2}));
end % if

%% Exact solution
if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
  problemData.pdExact = problemData;
end % if

%% Initial error computation.
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

fprintf('L2 errors of h, u1 w.r.t. the initial condition: %g, %g\n', ...
  computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{end, 1}, problemData.h0Cont, problemData.qOrd+1, problemData.basesOnQuad1D), ...
  computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 2}, problemData.u10Cont, problemData.qOrd+1, problemData.basesOnQuad2D));
end % function