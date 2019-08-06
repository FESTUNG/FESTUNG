% Compute integrals on the reference triangle, whose integrands 
% consist of all permutations of derivatives of two basis functions.

%===============================================================================
%> @file ./core/integrateRefElemDPhiDPhi.m
%>
%> @brief Compute integrals on the reference triangle, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions.
%>
%> It computes third order tensor as described in [KnabnerAngermann2003, p.79/80]
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @param  qOrd (optional) The order of the quadrature rule to be used.
%> @retval ret  The computed matrix @f$[N\times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Andreas Rupp, 2018.
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
function ret = integrateRefElemDPhiDPhi(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); end
[~,~,W] = quadRule2D(qOrd);
ret = zeros(N,N,3); % [N x N]
for i = 1 : N
  for j = 1 : N
    ret(i, j, 1) = sum( W' .* basesOnQuad.gradPhi2D{qOrd}(:, i,1) .* basesOnQuad.gradPhi2D{qOrd}(:, j,1) );
    ret(i, j, 2) = sum( W' .* basesOnQuad.gradPhi2D{qOrd}(:, i,1) .* basesOnQuad.gradPhi2D{qOrd}(:, j,2) ) ...
                 + sum( W' .* basesOnQuad.gradPhi2D{qOrd}(:, i,2) .* basesOnQuad.gradPhi2D{qOrd}(:, j,1) );
    ret(i, j, 3) = sum( W' .* basesOnQuad.gradPhi2D{qOrd}(:, i,2) .* basesOnQuad.gradPhi2D{qOrd}(:, j,2) );
  end % for
end % for
end % function
