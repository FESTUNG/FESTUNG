% Evaluates all permutations of two one-dimensional basis functions and the
% spatial derivative of a one-dimensional basis function in all quadrature 
% points on the unit interval.

%===============================================================================
%> @file integrateRefElem1DDphiPhiPhiPerQuad.m
%>
%> @brief Evaluates all permutations of two one-dimensional basis functions and
%>        the spatial derivative of a one-dimensional basis function in all 
%>        quadrature points on the unit interval.
%===============================================================================
%>
%> @brief Evaluates all permutations of two one-dimensional basis functions and
%>        the spatial derivative of a one-dimensional basis function in all 
%>        quadrature points on the unit interval.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{G}} 
%>    \in \mathbb{R}^{N\times N\times N\times R}@f$, which is defined by
%> @f[
%> \left[\hat{\mathsf{G}}\right]_{i,j,l,r} \;:=\;
%>   w_r \, \partial_{\hat{x}} \hat{\phi}_i(q_r) \,
%>   \hat{\phi}_j(q_r) \, \hat{\phi}_l(q_r)
%> @f]
%> with the quadrature points @f$q_r@f$ and weights @f$w_r@f$ given by
%> <code>quadRule1D()</code>
%>
%> @param  N            The local number of degrees of freedom
%> @param  qOrd         The order of the quadrature rule to be used
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D and gradPhi1D.
%> @retval ret  The computed array @f$[N\times N\times N\times R]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function ret = integrateRefElem1DDphiPhiPhiPerQuad(N, qOrd, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
[Q, W] = quadRule1D(qOrd); R = length(W);
ret = { zeros(N, N, N, R), zeros(N, N, N, R), zeros(N, N, N, R) }; % [N x N x N x R]
for i = 1 : N
  for j = 1 : N
    for l = 1 : N
      ret{1}(i, j, l, :) = W.' .* ( basesOnQuad.gradPhi1D{qOrd}(:, i) .* basesOnQuad.phi1D{qOrd}(:, j) .* basesOnQuad.phi1D{qOrd}(:, l) );
      ret{2}(i, j, l, :) = (W .* Q).' .* ( basesOnQuad.gradPhi1D{qOrd}(:, i) .* basesOnQuad.phi1D{qOrd}(:, j) .* basesOnQuad.phi1D{qOrd}(:, l) );
    end % for l
  end % for j
end % for i
end % function