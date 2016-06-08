% Fills the problem's data structures with initial data (if applicable).

%===============================================================================
%> @file template/initializeProblem.m
%>
%> @brief Fills the problem's data structures with initial data (if applicable).
%===============================================================================
%>
%> @brief Fills the problem's data structures with initial data (if applicable).
%>
%> This routine is called after template/preprocessProblem.m.
%>
%> Before entering the main loop of the solution algorithm, this routine
%> fills the problem's data structures with initial data.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by configureProblem() and 
%>                      preprocessProblem(). @f$[\text{struct}]@f$
%>
%> @retval problemData  The input struct enriched with initial data.
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
function pd = initializeProblem(pd)
p = pd.p;
K = pd.K;
N = pd.N;

%% Primary unknowns H, uH, vH (each K*N).
if pd.isSolutionAvail
  h0Cont = @(x1,x2) pd.xiCont(x1,x2,pd.t0) - pd.zbCont(x1,x2);
  uH0Cont = @(x1,x2) h0Cont(x1,x2) .* pd.uCont(x1,x2,pd.t0);
  vH0Cont = @(x1,x2) h0Cont(x1,x2) .* pd.vCont(x1,x2,pd.t0);
else
  h0Cont = @(x1,x2)  -pd.zbCont(x1,x2);
  uH0Cont = @(x1,x2) zeros(size(x1));
  vH0Cont = @(x1,x2) zeros(size(x1));
end % if
  
pd.cDisc = zeros(K,N,3);
pd.cDisc(:,:,1) = projectFuncCont2DataDisc(pd.g, h0Cont, 2*p, pd.refElemPhiPhi, pd.basesOnQuad);
pd.cDisc(:,:,2) = projectFuncCont2DataDisc(pd.g, uH0Cont, 2*p, pd.refElemPhiPhi, pd.basesOnQuad);
pd.cDisc(:,:,3) = projectFuncCont2DataDisc(pd.g, vH0Cont, 2*p, pd.refElemPhiPhi, pd.basesOnQuad);

for i = 1 : length(pd.slopeLimList)
  switch pd.slopeLimList{i}
    case 'H'
      dataV0T = pd.ramp(pd.t0 / 86400) * pd.g.markV0TbdrRI .* xiRivV0T + (pd.g.markV0TbdrOS & ~pd.g.markV0TbdrRI) .* xiOSV0T;
      pd.cDisc(:,:,1) = applySlopeLimiterDisc(pd.g, pd.cDisc(:,:,1), pd.g.markV0TbdrRI | pd.g.markV0TbdrOS, dataV0T, pd.globM, pd.globMDiscTaylor, pd.typeSlopeLim);
    otherwise
      error('Unknown variable for slope limiting')
  end % switch
end % for

pd.cDisc(:,:,1) = correctMinValueExceedanceDisc(pd.cDisc(:,:,1), pd.sysMinValueCorrection, 0, pd.minValueHeight, 1000);

%% Visualize initial solution.
visualizeSolution(pd, 0);

%% Initialize waitbar.
if pd.isWaitbar
  str  = strcat( ['% done. Simulating refinement level =', ' ', num2str(pd.refinement), ', p =', ' ', num2str(p), '.' ]);
  pd.waitbar = waitbar(0, strcat([ 'Time stepping:', ' ', num2str(0), str ]));
end % if

%% Initialize time stepping
pd.isFinished = false;
pd.t = pd.t0;
end % function