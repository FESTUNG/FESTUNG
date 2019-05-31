% Performs all post-processing steps, such as error estimates of the final
% solution, etc.

%===============================================================================
%> @file ./cahn-hilliard/configureProblem.m
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%===============================================================================
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%>
%> This routine is called after the main loop.
%>
%> It can include things such as error estimates, an output operation of
%> the final solution, etc.
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

if strcmp(problemData.scenario,'conv-test')
  K = problemData.K;
  N = problemData.N;
  p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);
  
  sol = @(x,y) problemData.cCont(problemData.tEnd,x,y);
  cDisc = reshape(problemData.sysY(1:K*N), N, K)';
  problemData.L2Error = computeL2Error(problemData.g, cDisc, sol, qOrd, problemData.basesOnQuad);
  problemData.L2Error
end

if isfield(problemData,'maxResidual')
  fprintf('Newton failed %d time(s). Maximum residual was %.3e.\n', problemData.newtonFails, problemData.maxResidual);
end

if problemData.gong == true
  load gong.mat;
  sound(y);
end

end

