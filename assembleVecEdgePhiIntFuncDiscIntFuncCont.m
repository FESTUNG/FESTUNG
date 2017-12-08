% Assembles a vector containing integrals over edges of products of a basis 
% function with a discontinuous coefficient function and a continuous function.

%===============================================================================
%> @file assembleVecEdgePhiIntFuncDiscIntFuncCont.m
%>
%> @brief Assembles a vector containing integrals over edges of products of a 
%>        basis function with a discontinuous coefficient function and a 
%>        continuous function.
%===============================================================================
%>
%> @brief Assembles a vector @f$\mathbf{K}_\mathrm{N}@f$ containing integrals 
%>        over edges of products of a basis function with a discontinuous 
%>        coefficient function and a continuous function
%>        @f$g_\mathrm{N}(t, \mathbf{x})@f$.
%>
%> The vector @f$\mathbf{K}_\mathrm{N} \in \mathbb{R}^{KN}@f$ is defined
%> component-wise by
%> @f[
%> [\mathbf{K}_\mathrm{N}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} \sum_{l=1}^N D_{kl}(t)
%>  \int_{E_{kn}} \varphi_{ki} \varphi_{kl} g_\mathrm{N}(t) \mathrm{d}s\,.
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
%>   \int_{E_{kn}} \varphi_{ki} \varphi_{kl} g_\mathrm{N}(t) \mathrm{d}s =
%>   \frac{|E_{kn}|}{|\hat{E}_n|} \int_{\hat{E}_n}\hat{\varphi}_i\hat{\varphi}_l
%>   g_\mathrm{N}(t,\mathbf{F}_k(\hat{\mathbf{x}}))\mathrm{d}\hat{\mathbf{x}}\,.
%> @f]
%> Further transformation to the unit interval @f$[0,1]@f$ using the mapping 
%> @f$\hat{\mathbf{\gamma}}_n(s)@f$ as provided by <code>gammaMap()</code>
%> gives the component-wise formulation
%> @f[
%>  [\mathbf{K}_\mathrm{N}]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} |E_{kn}| \sum_{l=1}^N
%>  D_{kl}(t) \int_0^1 \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(s)
%>     \hat{\varphi}_{l} \circ \hat{\mathbf{\gamma}}_n(s)
%>     g_\mathrm{N}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(s))
%>     \mathrm{d}s \,.
%> @f]
%> This integral is then approximated using a 1D quadrature rule provided by
%> <code>quadRule1D()</code> 
%> @f[
%>  [\mathbf{K}_\mathrm{N}]_{(k-1)N+i} \approx
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} |E_{kn}| \sum_{l=1}^N
%>  D_{kl}(t) \sum_{r=1}^R \omega_r
%>     \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     \hat{\varphi}_{l} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     g_\mathrm{N}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(q_r)) \,,
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
%> @param dataDisc    A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @param areaE0Tbdr  (optional) argument to provide precomputed values
%>                    for the products of <code>markE0Tbdr</code>,
%>                    and <code>g.areaE0T</code>,
%>                    @f$[3 \text{ cell}]@f$
%> @retval ret        The assembled vector @f$[KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> Modified by Hennes Hajduk, 2016-04-06
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
function ret = assembleVecEdgePhiIntFuncDiscIntFuncCont(g, markE0Tbdr, dataDisc, funcCont, basesOnQuad, areaE0Tbdr)
% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [g.numT 3]}, mfilename, 'markE0Tbdr');
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(dataDisc, {'numeric'}, {'size', [g.numT NaN]}, mfilename, 'dataDisc');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

if nargin > 5
  ret = assembleVecEdgePhiIntFuncDiscIntFuncCont_WithAreaE0Tbdr(g, dataDisc, funcCont, basesOnQuad, areaE0Tbdr);
else
  ret = assembleVecEdgePhiIntFuncDiscIntFuncCont_NoAreaE0Tbdr(g, markE0Tbdr, dataDisc, funcCont, basesOnQuad);
end % if
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleVedEdgePhiIntFuncDiscIntFuncCont()
%> was called with a precomputed field areaE0Tbdr.
%
function ret = assembleVecEdgePhiIntFuncDiscIntFuncCont_WithAreaE0Tbdr(g, dataDisc, funcCont, basesOnQuad, areaE0Tbdr)
% Determine quadrature rule and mapping to physical element
[K, N] = size(dataDisc);  p = (sqrt(8*N+1)-3)/2;
qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd);

% Assemble vector
ret = zeros(K, N);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  funcAtQ = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  for i = 1 : N
    for l = 1 : N
      integral = funcAtQ * ( W .* basesOnQuad.phi1D{qOrd}(:,i,n)' .* basesOnQuad.phi1D{qOrd}(:,l,n)' )';
      ret(:,i) = ret(:,i) + areaE0Tbdr{n} .* dataDisc(:,l) .* integral;
    end % for
  end % for
end % for
ret = reshape(ret',K*N,1);
end % function
%
%===============================================================================
%> @brief Helper function for the case that assembleVedEdgePhiIntFuncDiscIntFuncCont()
%> was called with no precomputed field areaE0Tbdr.
%
function ret = assembleVecEdgePhiIntFuncDiscIntFuncCont_NoAreaE0Tbdr(g, markE0Tbdr, dataDisc, funcCont, basesOnQuad)
% Determine quadrature rule and mapping to physical element
[K, N] = size(dataDisc);  p = (sqrt(8*N+1)-3)/2;
qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd);

% Assemble vector
ret = zeros(K, N);
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  funcAtQ = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  Kkn = markE0Tbdr(:,n) .* g.areaE0T(:,n);
  for i = 1 : N
    for l = 1 : N
      integral = funcAtQ * ( W .* basesOnQuad.phi1D{qOrd}(:,i,n)' .* basesOnQuad.phi1D{qOrd}(:,l,n)' )';
      ret(:,i) = ret(:,i) + Kkn .* dataDisc(:,l) .* integral;
    end % for
  end % for
end % for
ret = reshape(ret',K*N,1);
end % function
