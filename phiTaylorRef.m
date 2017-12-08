% Evaluates the i-th Taylor basis function in on all elements in physical domain
% according to a given point in the reference element.
%
%===============================================================================
%> @file phiTaylorRef.m
%>
%> @brief Evaluates the i-th Taylor basis function in on all elements in
%>        physical domain according to a given point in the reference element.
%===============================================================================
%>
%> @brief Evaluates the i-th Taylor basis function in on all elements in
%>        physical domain according to a given point in the reference element.
%>
%> It computes 
%> @f$\hat{\phi}(\hat{\mathbf{x}}) = \phi \circ \mathbf{F}_{1:K}(\hat{\mathbf{x}})@f$
%> using <code>phiTaylorPhy()</code>.
%> 
%> @param  g      The lists describing the geometric and topological 
%>                properties of a triangulation (see 
%>                <code>generateGridData()</code>) 
%>                @f$[1 \times 1 \text{ struct}]@f$
%> @param  i      The index of the basis function.
%> @param  hatX1  A list of @f$\hat{x}^1@f$ coordinates.
%>                @f$[1 \times n_\mathrm{Points}@f$]
%> @param  hatX2  A list of @f$\hat{x}^2@f$ coordinates.
%>                @f$[1 \times n_\mathrm{Points}@f$]
%> @retval ret    The @f$i@f$-th basis function in all points specified by
%>                <code>hatX1</code>, <code>hatX2</code>.
%>                It holds <code>size(hatX1) == size(hatX2) == size(ret)</code>.
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
function ret = phiTaylorRef(g, i, hatX1, hatX2)
ret = phiTaylorPhy(g, i, g.mapRef2Phy(1, hatX1, hatX2), g.mapRef2Phy(2, hatX1, hatX2));
end

