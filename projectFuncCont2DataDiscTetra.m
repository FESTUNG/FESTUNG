% Compute the DG/modal basis representation of an algebraic function.

%===============================================================================
%> @file projectFuncCont2DataDiscTetra.m
%>
%> @brief Compute the DG/modal basis representation of an algebraic function.
%===============================================================================
%>
%> @brief Compute the DG/modal basis representation of an algebraic function by
%>        performing the @f$L^2@f$-projection.
%>
%> The DG/modal basis representation of a algebraic function @f$d(t)@f$,
%> is given as
%> @f[
%>   d_h(t, \mathbf{x}) = \sum_{j=1}^N D_{kj}(t) \varphi_{kj}(\mathbf{x}) \,,
%> @f]
%> such that @f$d_h(t) \in \mathbb{P}_d (\mathcal{T}_h)@f$ is an adequate
%> approximation of @f$d(t)@f$.
%>
%> To produce @f$d_h(t)@f$ we use the @f$L^2@f$-projection defined locally for
%> @f$T_k \in \mathcal{T}_h@f$ by
%> @f[
%>  \forall w_h \in \mathbb{P}_d(T_k), \quad
%>  \int_{T_k} w_h d_h(t) = \int_{T_k} w_h d(t) \,.
%> @f]
%> Choosing @f$w_h = \varphi_{ki} \text{ for } i \in \{1,...,N\}@f$ and using
%> a mapping 
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$ from 
%> the reference element @f$\hat{T}@f$ we obtain
%> @f[
%>  \sum_{j=1}^N D_{kj}(t) \int_{\hat{T}} \hat{\varphi}_i(\hat{\mathbf{x}})
%>  \hat{\varphi}_j(\hat{\mathbf{x}}) \mathrm{d} \hat{\mathbf{x}} =
%> \int_{\hat{T}} \hat{\varphi}_i(\hat{\mathbf{x}}) 
%> d(t, \mathbf{F}_k(\hat{\mathbf{x}})) \mathrm{d} \hat{\mathbf{x}} \,.
%> @f]
%> Written in matrix form, this is equivalent to
%> @f[
%>  \hat{\mathsf{M}} \begin{bmatrix} D_{k1} \\ \vdots \\ D_{kN} \end{bmatrix} =
%>  \int_{\hat{T}} \begin{bmatrix} 
%>    \hat{\varphi}_1(\hat{\mathbf{x}}) d(t, \mathbf{F}_k(\hat{\mathbf{x}})) \\
%>    \vdots \\
%>    \hat{\varphi}_N(\hat{\mathbf{x}}) d(t, \mathbf{F}_k(\hat{\mathbf{x}}))
%>  \end{bmatrix} \mathrm{d} \hat{\mathbf{x}} \,,
%> @f]
%> with local mass matrix on the reference element 
%> @f$\hat{M} \in \mathbb{R}^{N\times N}@f$ as defined in
%> <code>integrateRefElemPhiPhi</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  N          Number of local degrees of freedom
%> @param  qOrd       The order of the quadrature rule.
%> @param  globM      Global mass matrix @f$\mathsf{M}@f$ as provided
%>                    by <code>assembleMatElemPhiPhi()</code>.
%>                    @f$[KN \times KN]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval The representation matrix of the DG/modal basis representation.
%>         @f$[K \times N]@f$
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
function dataDisc = projectFuncCont2DataDiscTetra(g, funcCont, N, qOrd, globM, basesOnQuad)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(globM, {'numeric'}, {'size', [g.numT * N, g.numT * N]}, mfilename, 'globM');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad');
[Q, W] = quadRule1D(max(qOrd,1)); [Q1, Q2] = meshgrid(Q); W = W' * W;
Q1 = Q1(:)'; Q2 = Q2(:)'; W = W(:);
rhs = (funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2)) .* ...
        (g.detJ0T{1} * ones(size(Q1)) + g.detJ0T{2} * Q1 + g.detJ0T{3} * Q2)) * ...
      (repmat(W, 1, N) .* basesOnQuad.phi2D{qOrd}(:,1:N));
dataDisc = reshape(globM \ reshape(rhs.', g.numT * N, 1), N, g.numT).';
end % function
