% Assembles two matrices, each containing integrals of products of a 
% one-dimensional basis  function with a (spatial) derivative of a 
% two-dimensional basis function.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices, each containing integrals of products of a 
%>        one-dimensional basis function with a (spatial) derivative of a 
%>        two-dimensional basis function.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{H}^m, m \in \{1,2\}@f$
%>        containing integrals of products of a one-dimensional basis function 
%>        with a (spatial) derivative of a two-dimensional basis function.
%>
%> The matrices @f$\mathsf{H}^m \in 
%>   \mathbb{R}^{KN\times \overline{K}\overline{N}}@f$ are block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{H}^m]_{(k-1)N+i,(\overline{k}-1)\overline{N}+j} = \int_{T_k} 
%>      \partial_{x^m}\varphi_{ki} \phi_{\overline{k}j} \mathrm{d}\mathbf{x}\,.
%> @f]
%> All other entries are zero.
%> For the implementation, the element integrals are backtransformed to the
%> reference square @f$\hat{T} = [0,1]^2@f$ using a mapping
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$.
%>
%> @param  g2D        The lists describing the geometric and topological 
%>                    properties of a two-dimensional triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a matching one-dimensional triangulation 
%>                    (see <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhi1D Local matrix @f$\hat{\mathsf{H}}@f$ as provided
%>                    by <code>integrateRefElemQuadriDphiPhi1D()</code>.
%>                    @f$[N \times \overline{N} \times 2]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
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
function ret = assembleMatElemQuadriDphiPhi1D(g2D, g1D, refElemDphiPhi1D)
K = g2D.numT; barK = g1D.numT; N = size(refElemDphiPhi1D{1});
ret = { sparse(K*N(1), barK*N(2)), sparse(K*N(1), barK*N(2)) };
for m = 1 : 2
  for s = 1 : 3
    ret{m} = ret{m} + ...
             kron(bsxfun(@times, g1D.markT2DT, g2D.J0T{s}(:,3-m,3-m)), refElemDphiPhi1D{s}(:,:,  m)) - ...
             kron(bsxfun(@times, g1D.markT2DT, g2D.J0T{s}(:,3-m,  m)), refElemDphiPhi1D{s}(:,:,3-m));
  end % for s
end % for m
end % function