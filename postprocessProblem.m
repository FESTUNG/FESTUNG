% Performs all post-processing steps. Error evaluation for analytical problems.

%===============================================================================
%> @file sweVert/postprocessProblem.m
%>
%> @brief Performs all post-processing tasks. Error evaluation for analytical problems.
%===============================================================================
%>
%> @brief Performs all post-processing tasks. Error evaluation for analytical problems.
%>
%> This routine is called after the main loop.
%>
%> @param  problemData  A struct with problem parameters and solution
%>                      vectors. @f$[\text{struct}]@f$
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
function problemData = postprocessProblem(problemData)
if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
  htEndCont = @(x1) problemData.hCont(problemData.tEnd, x1);
  u1tEndCont = @(x1,x2) problemData.u1Cont(problemData.tEnd, x1, x2);
  u2tEndCont = @(x1,x2) problemData.u2Cont(problemData.tEnd, x1, x2);

  problemData.error = [ computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{1, 1}, ...
                            htEndCont, problemData.qOrd + 1, problemData.basesOnQuad1D), ...
                        execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDiscRK{1, 2}, ...
                            u1tEndCont, problemData.qOrd + 1, problemData.basesOnQuad2D), ...
                        execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDiscRK{1, 3}, ...
                            u2tEndCont, problemData.qOrd + 1, problemData.basesOnQuad2D) ];

  fprintf('L2 errors of cDisc w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
end % if
if problemData.isVisGrid, execin('darcyVert/visualizeGridTrap', problemData.g); end
end % function

