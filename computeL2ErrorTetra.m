% Compute the L2-error of a DG/modal basis representation on quadrilaterals 
% for a given analytical solution.

%===============================================================================
%> @file computeL2ErrorTetra.m
%>
%> @brief Compute the L2-error of a DG/modal basis representation on
%>        quadrilaterals for a given analytical solution.
%===============================================================================
%>
%> @brief Compute the L2-error of a DG/modal basis representation on
%>        quadrilaterals for a given analytical solution.
%>
%> The discretization error @f$\|c_h(t)-c(t)\|_{L^2(\Omega)}@f$ at time @f$t@f$
%> gives the @f$L^2@f$-norm of the difference between the discrete solution 
%> @f$c_h(t)@f$ and the analytical solution @f$c(t)@f$.
%>
%> To compute the discretization error, we must evaluate
%> @f[
%>  \|c_h(t)-c(t)\|^2_{L^2(\Omega)} = 
%>    \sum_{T_k\in\mathcal{T}_h} \int_{T_k} (c_h(t)-c(t))^2 \mathrm{d}\mathbf{x}.
%> @f]
%> where the DG/modal basis representation of the discrete solution @f$c(t)@f$,
%> is given as
%> @f[
%>   c_h(t, \mathbf{x}) = \sum_{j=1}^N C_{kj}(t) \varphi_{kj}(\mathbf{x}) \,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g.,
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc   The representation matrix of the DG/modal basis
%>                    representation @f$[K \times N]@f$
%> @param  funcCont   A function handle for the continuous solution
%> @param  qOrd       The order of the quadrature rule to be used.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret The discretization error @f$[\text{scalar}]@f$
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
function err = computeL2ErrorTetra(g, dataDisc, funcCont, qOrd, basesOnQuad)
% Determine quadrature rule and physical coordinates
[Q, W] = quadRule1D(qOrd); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:)'; Q2 = Q2(:)'; W = W(:)';
X1 = g.mapRef2Phy(1,Q1,Q2); X2 = g.mapRef2Phy(2,Q1,Q2);
N = size(dataDisc,2);

% Evaluate analytical and discrete function
fExOnQuadPts = funcCont(X1, X2);
fApprxOnQuadPts = dataDisc * basesOnQuad.phi2D{qOrd}(:,1:N).'; % [K x R] = [K x N] * [N x R]

% Compute error
err = sqrt(dot((fApprxOnQuadPts - fExOnQuadPts).^2 * W.', g.areaT));
end  % function

