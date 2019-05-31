% Compute the DG/modal basis representation of a function evaluated in the
% physical quadrature points of a grid.

%===============================================================================
%> @file
%>
%> @brief Compute the DG/modal basis representation of a function evaluated 
%>				in the physical quadrature points of a grid.
%===============================================================================
%>
%> @brief Compute the DG/modal basis representation of a function evaluated 
%>				in the physical quadrature points of a grid by
%>        performing the @f$L^2@f$-projection.
%>
%> The DG/modal basis representation of a function @f$d(t,\mathbf{x})@f$
%> is given as
%> @f[
%>   d_h(t, \mathbf{x}) = \sum_{j=1}^N D_{kj}(t) \varphi_{kj}(\mathbf{x}) \,,
%> @f]
%> such that @f$d_h(t,\mathbf{x}) \in \mathbb{P}_d (\mathcal{T}_h)@f$ is an 
%> adequate approximation of @f$d(t,\mathbf{x})@f$.\n
%> In this case we do not know an analytical espression for @f$d(t)@f$, but its
%> values in appropriate quadrature points.
%>
%> To produce @f$d_h(t,\mathbf{x})@f$ we use the @f$L^2@f$-projection defined
%> locally for @f$T_k \in \mathcal{T}_h@f$ by
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
%> <code>integrateRefElemPhiPhi</code>.\n
%> Since the right hand side is approximated by quadrature, we have all we need.
%>
%> @param  dataQ0T    A coefficient matrix representing the function @f$d(t)@f$
%>										where each row contains the evaluations of @f$d(t)@f$
%>										in the quadrature points of the corresponding element
%>										@f$[K \times N]@f$
%> @param  ord        The order of the quadrature rule provided by 
%>                    <code>quadRule2D()</code>
%> @param refElemPhiPhi Local matrix @f$\hat{\mathsf{M}}@f$ as provided
%>                    by <code>integrateRefElemPhiPhi()</code>.
%>                    @f$[N \times N]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%>
%> @retval The representation matrix of the DG/modal basis representation.
%>         @f$[K \times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function dataDisc = projectDataQ0T2DataDisc(dataQ0T, ord, refElemPhiPhi, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad');
ord = max(ord,1); [~, ~, W] = quadRule2D(ord);
N = size(refElemPhiPhi, 1);
rhs = dataQ0T * (repmat(W.', 1, N) .* basesOnQuad.phi2D{ord});
dataDisc = rhs / refElemPhiPhi;
end % function
