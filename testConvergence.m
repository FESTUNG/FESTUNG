% Driver routine used to trigger testConvergence() implementations in problem 
% folders.

%===============================================================================
%> @file
%>
%> @brief Driver routine used to trigger testConvergence() implementations in 
%>        problem folders.
%===============================================================================
%>
%> @brief Driver routine used to trigger testConvergence() implementations in 
%>        problem folders.
%>
%> It solves a given testcase with continuously refined mesh size and/or
%> time step size and saves the obtained error estimates to compute
%> experimental orders of convergence.
%>
%> For that, it calls the implementation of testConvergence() in the problem
%> folder of the given problem name.
%>
%> @param  problemName  The name of the problem to be used.
%> @param  varargin     Arbitrary arguments that are passed to the called
%>                      implementation of testConvergence().
%>
%> @retval err          A cell array of vectors, with errors from the
%>                      computations on all space/time levels for each
%>                      polynomial degree.
%> @retval eoc          A cell array with the corresponding experimental
%>                      orders of convergence.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Balthasar Reuter, 2017.
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
function [err, eoc] = testConvergence(problemName, varargin)
oldpath = addpath([pwd filesep 'core']);
fn_testConvergence = getFunctionHandle([problemName filesep 'doConvergenceTest']);
path(oldpath);
[err, eoc] = fn_testConvergence(varargin{:});
end % function
