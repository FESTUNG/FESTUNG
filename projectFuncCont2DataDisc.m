% Compute the DG/modal basis representation of an algebraic function.
%
%===============================================================================
%> @file projectFuncCont2DataDisc.m
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
%> an affine mapping 
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$ from 
%> the reference triangle @f$\hat{T} = \{(0,0), (1,0), (0,1)\}@f$, defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,,
%> @f]
%> we obtain
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
%> with local mass matrix on the reference triangle 
%> @f$\hat{M} \in \mathbb{R}^{N\times N}@f$ as defined in
%> <code>integrateRefElemPhiPhi</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  ord        The order of the quadrature rule provided by 
%>                    <code>quadRule2D()</code>
%> @param refElemPhiPhi Local matrix @f$\hat{\mathsf{M}}@f$ as provided
%>                    by <code>integrateRefElemPhiPhi()</code>.
%>                    @f$[N \times N]@f$
%> @retval The representation matrix of the DG/modal basis representation.
%>         @f$[K \times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function dataDisc = projectFuncCont2DataDisc(g, funcCont, ord, refElemPhiPhi)
global gPhi2D
ord = max(ord,1);  [Q1, Q2, W] = quadRule2D(ord);
K = g.numT; N = size(refElemPhiPhi, 1);
F1 = @(X1, X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
F2 = @(X1, X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
rhs = funcCont(F1(Q1, Q2), F2(Q1, Q2)) * (repmat(W.', 1, N) .* gPhi2D{ord});
dataDisc = rhs / refElemPhiPhi;
end % function
