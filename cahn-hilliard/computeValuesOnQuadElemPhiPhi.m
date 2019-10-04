% Assembles a third prder tensor containing the values of the product of
% two basis functions at the quardature points

%===============================================================================
%> @file assembleMatElemDphiPhiFuncDisc.m
%>
%> @brief Assembles a third prder tensor containing the values of the product
%         of two basis functions at the quardature points.
%===============================================================================
%>
%> @brief Assembles a third prder tensor containing the values of the product
%         of two basis functions at the quardature points.
%>
%> The tensor @f$\mathsf{VoQ} \in \mathbb{R}^{N \times N \times Q}@f$ are
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{VoQ}^m]_{i,j,r} = \varphi_i(x_r) \varphi_j(x_r)\,.
%> @f]
%> Here, q is the number of quadrature point and x_r is the r-th quadrature
%> point.
%>
%> @param  N          The number of local degrees of freedom. For polynomial
%>                    order @f$p@f$, it is given as @f$N = (p+1)(p+2)/2@f$
%>                    @f$[\text{scalar}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret        The assembled matrix
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = computeValuesOnQuadElemPhiPhi(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2; qOrd = max(2*p, 1); end
[~,~,W] = quadRule2D(qOrd);
R = length(W);

% Compute value of phi_i * phi_j on quadrature point r
ret = zeros(size(basesOnQuad.phi2D{qOrd}, 2), size(basesOnQuad.phi2D{qOrd}, 2), size(basesOnQuad.phi2D{qOrd}, 1));
for r = 1 : R
  ret(:,:,r) = W(r) * kron( basesOnQuad.phi2D{qOrd}(r,:), basesOnQuad.phi2D{qOrd}(r,:).' );
end % for r

end % function computeValuesOnQuadElemPhiPhi