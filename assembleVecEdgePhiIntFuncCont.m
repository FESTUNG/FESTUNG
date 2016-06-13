% Assembles a vector containing integrals over edges of products of a basis 
% function with a continuous function.

%===============================================================================
%> @file assembleVecEdgePhiIntFuncCont.m
%>
%> @brief Assembles a vector containing integrals over edges of products of a 
%>        basis function with a continuous function.
%===============================================================================
%>
%> @brief Assembles a vector @f$\mathbf{K}_\mathrm{D}@f$ containing integrals 
%>        over edges of products of a basis function with a continuous function
%>        @f$c_\mathrm{D}(t, \mathbf{x})@f$.
%>
%> The vector @f$\mathbf{K}_\mathrm{D} \in \mathbb{R}^{KN}@f$ is defined
%> component-wise by
%> @f[
%> [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \frac{1}{|E_{kn}|} \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s\,.
%> @f]
%> For the implementation, the integrals are backtransformed to the
%> reference triangle @f$\hat{T} = \{(0,0), (1,0), (0,1)\}@f$ using an affine
%> mapping @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$
%> defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,.
%> @f]
%> This allows to reformulate 
%> @f[
%>   \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s =
%>   \frac{|E_{kn}|}{|\hat{E}_n|} \int_{\hat{E}_n} \hat{\varphi}_i 
%>   c_\mathrm{D}(t,\mathbf{F}_k(\hat{\mathbf{x}}))\mathrm{d}\hat{\mathbf{x}}\,.
%> @f]
%> Further transformation to the unit interval @f$[0,1]@f$ using the mapping 
%> @f$\hat{\mathbf{\gamma}}_n(s)@f$ as provided by <code>gammaMap()</code>
%> gives the component-wise formulation
%> @f[
%>  [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \int_0^1 \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(s)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(s))
%>     \mathrm{d}s \,.
%> @f]
%> This integral is then approximated using a 1D quadrature rule provided by
%> <code>quadRule1D()</code> 
%> @f[
%>  [\mathbf{K}_\mathrm{D}]_{(k-1)N+i} \approx
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \sum_{r=1}^R \omega_r \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(q_r)) \,,
%> @f]
%> allowing to vectorize over all triangles.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0Tbdr <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the vector entries should be
%>                    assembled @f$[K \times 3]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  N          The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret        The assembled vector @f$[KN]@f$
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
function ret = assembleVecEdgePhiIntFuncCont(g, markE0Tbdr, funcCont, N, basesOnQuad)
% Determine quadrature rule and mapping to physical element
p = (sqrt(8*N+1)-3)/2;
qOrd = 2*p+1; [Q, W] = quadRule1D(qOrd);
Q2X1   = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
Q2X2   = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));

% Check function arguments that are directly used
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(markE0Tbdr, {'logical', 'numeric'}, {'size', [g.numT 3]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

% Assemble vector
ret = zeros(g.numT, N);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  cDn = funcCont(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
  for i = 1 : N
    ret(:,i) = ret(:,i) + markE0Tbdr(:,n) .* ( cDn * (W' .* basesOnQuad.phi1D{qOrd}(:,i,n)) );
  end % for
end % for
ret = reshape(ret', g.numT*N, 1);
end % function
