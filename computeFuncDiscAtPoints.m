% Assembles a matrix containing the values of a discrete function given by 
% a representation matrix evaluated in all points, for which the evaluated
% basis functions are given.
%
%===============================================================================
%> @file computeFuncDiscAtPoints.m
%>
%> @brief Assembles a matrix containing the values of a discrete function
%>        given by a representation matrix evaluated in all points, for
%>        which the evaluated basis functions are given.
%===============================================================================
%>
%> @brief Assembles a matrix containing the values of a discrete function
%>        given by a representation matrix evaluated in all points, for
%>        which the evaluated basis functions are given.
%>
%> @param  funcDisc    The representation matrix of the discrete function.
%>                     @f$[K \times N]@f$
%> @param  phiAtPoints The evaluated basis functions in each point of each
%>                     triangle. @f$[K\times n_\mathrm{Points}\times N]@f$
%> @retval ret         The assembled matrix @f$[K \times n_\mathrm{Points}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = computeFuncDiscAtPoints(funcDisc, phiAtPoints)
% Extract and verify dimensions
nPoints = size(phiAtPoints, 2);
[K, N] = size(funcDisc);
assert(K == size(phiAtPoints, 1), 'Number of elements does not match!');
assert(size(phiAtPoints, 3) == N, 'Number of DOF does not match!');

% Evaluate function in all given points
ret = zeros(K, nPoints);
for i = 1 : nPoints
  ret(:, i) = sum(funcDisc .* squeeze(phiAtPoints(:, i, :)), 2);
end %for

end % function
