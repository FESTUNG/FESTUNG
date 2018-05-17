% Fills the problem's data structures with initial data.

%===============================================================================
%> @file
%>
%> @brief Fills the problem's data structures with initial data.
%===============================================================================
%>
%> @brief Fills the problem's data structures with initial data.
%>
%> This routine is called after darcy_2dv/preprocessProblem.m.
%>
%> Before entering the main loop of the solution algorithm, this routine
%> fills the problem's data structures with initial data.
%>
%> It projects the initial condition and initializes the solution vector.
%> If <tt>problemData.isHotstart</tt> is set to true, it loads the current
%> state from the file specified in <tt>problemData.hotstart</tt>.
%>
%> If visualization is enabled (i.e., <tt>problemData.isVisSol == true</tt>)
%> the initial state is written to a visualization file.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by darcy_2dv/configureProblem.m and 
%>                      darcy_2dv/preprocessProblem.m. @f$[\text{struct}]@f$
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
K = problemData.g.numT;
N = problemData.N;

%% Initial data
if problemData.isHotstart
  hotstart = load(problemData.hotstartFile);
  assert(all(isfield(hotstart, {'t', 'hDisc', 'q1Disc', 'q2Disc'})), 'Hotstart file must contain t, hDisc, q1Disc, q2Disc')
  assert(all(problemData.g.numT == [size(hotstart.hDisc, 1), size(hotstart.q1Disc, 1), size(hotstart.q2Disc, 1)]), ...
    'Number of elements in hotstart file does not match!')
  fprintf('Loaded hotstart data from "%s" at time level t=%g.\n', problemData.hotstartFile, hotstart.t);
  
  hDisc = zeros(K, N);
  q1Disc = zeros(K, N);
  q2Disc = zeros(K, N);
  
  hotstartN = min(N, size(hotstart.hDisc, 2));
  hDisc(:, 1:hotstartN) = hotstart.hDisc(:, 1:hotstartN);
  q1Disc(:, 1:hotstartN) = hotstart.q1Disc(:, 1:hotstartN);
  q2Disc(:, 1:hotstartN) = hotstart.q2Disc(:, 1:hotstartN);
else
  if isfield(problemData, 'hCont')
    h0Cont = @(x1,x2) problemData.hCont(0,x1,x2);
  elseif isfield(problemData, 'h0Cont')
    h0Cont = problemData.h0Cont;
  else
    warning('No initial data for h given. Initializing with zeros.');
    h0Cont = @(x1,x2) zeros(size(x1));
  end % if

  if isfield(problemData, 'q1Cont')
    q10Cont = @(x1,x2) problemData.q1Cont(0,x1,x2);
  elseif isfield(problemData, 'q10Cont')
    q10Cont = problemData.q10Cont;
  else
    warning('No initial data for q1 given. Initializing with zeros.');
    q10Cont = @(x1,x2) zeros(size(x1));
  end % if

  if isfield(problemData, 'q2Cont')
    q20Cont = @(x1,x2) problemData.q2Cont(0,x1,x2);
  elseif isfield(problemData, 'q20Cont')
    q20Cont = problemData.q20Cont;
  else
    warning('No initial data for q2 given. Initializing with zeros.');
    q20Cont = @(x1,x2) zeros(size(x1));
  end % if

  hDisc = projectFuncCont2DataDiscTetra(problemData.g, h0Cont, problemData.qOrd, ...
                                        problemData.globM, problemData.basesOnQuad);
  q1Disc = projectFuncCont2DataDiscTetra(problemData.g, q10Cont, problemData.qOrd, ...
                                         problemData.globM, problemData.basesOnQuad);
  q2Disc = projectFuncCont2DataDiscTetra(problemData.g, q20Cont, problemData.qOrd, ...
                                         problemData.globM, problemData.basesOnQuad);
                                       
  fprintf('L2 error w.r.t. the initial condition: %g\n', ...
    computeL2ErrorTetra(problemData.g, hDisc, h0Cont, problemData.qOrd+1, problemData.basesOnQuad));
end % if

problemData.sysY = [ reshape(q1Disc', K * N, 1) ; ...
                     reshape(q2Disc', K * N, 1) ; ...
                     reshape(hDisc', K * N, 1) ];

%% Visualization of inital condition.
if problemData.isVisSol
  cLagr = { projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)'), ...
            projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(1 : K*N), N, K)'), ...
            projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(K*N+1 : 2*K*N), N, K)') };
  visualizeDataLagrTetra(problemData.g, cLagr, {'h', 'q1', 'q2'}, problemData.outputBasename, 0, problemData.outputTypes, struct('q', {{'q1','q2'}}));
end % if

end % function