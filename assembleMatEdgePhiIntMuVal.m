% Assembles a matrix containing integrals of products of flux function and test 
% function.

%===============================================================================
%> @file assembleMatEdgePhiIntMuVal.m
%>
%> @brief TODO
%===============================================================================
%>
%> @brief TODO
%>
%> 
%> 
%> @f[
%>   [\mathsf{S}]_{(k-1)N+i,(\bar{k}-1)\bar{N}+j} =\sum_{E_{kn} \in \partial{{T_{k}}} \cap \ensuremath{\mathcal{E}}_{\text{int}} }  \int_{E_{kn}} ({\boldsymbol{u}}\cdot \boldsymbol{\nu}_{kn}) \,\mu_{knj}\, \varphi_{ki}\,  \text{d}s\,.
%> @f]
%> 
%> 
%> 
%> @f[
%> \mathsf{S}_{E_{kn}} = \int_{E_{kn}} ({\boldsymbol{u}} \cdot \boldsymbol{\nu}_{kn})
%> \begin{bmatrix}
%> \mu_{kn1} \, \varphi_{k1} & \ldots & \mu_{kn\bar{N}}  \, \varphi_{k1} \\
%> \vdots           & \ddots & \vdots \\
%> \mu_{kn1}  \, \varphi_{kN} & \ldots & \mu_{kn\bar{N}}  \, \varphi_{kN} 
%> \end{bmatrix}
%> \text{d}s\,.
%> @f]
%> 
%> 
%> 
%> @f[
%>   [\mathsf{S}_{\text{out}}]_{(k-1)N+i,(\bar{k}-1)\bar{N}+j} 
%>   = \sum_{E_{kn}\in\partial{T_{k}}\cap\ensuremath{\mathcal{E}}_{\text{out}}} \int_{E_{kn}} ({\boldsymbol{u}}\cdot\boldsymbol{\nu}_{kn}) \, \mu_{knj} \, \varphi_{ki} \,\text{d}s\,.
%> @f]
%> 
%> 

%> @f{align*}
%> \int_{E_{kn}} ({\boldsymbol{u}} \cdot \boldsymbol{\nu}_{kn}) \, \varphi_{ki} \, \mu_{knj} \, \text{d}s 
%> &= \ensuremath{|E_{kn}|} \int_0^1 \left(({\boldsymbol{u}}(t) \circ \boldsymbol{F}_{k} \circ \boldsymbol{\hat{\gamma}}_{n}(s)) \cdot \boldsymbol{\nu}_{kn} \right) \, \hat{\varphi}_{i} \circ \boldsymbol{\hat{\gamma}}_{n}(s) \, \hat{\mu}_{j} \circ \hat{\beta}_{kn}(s) \, \text{d}s \\
%> &\approx \sum_{r=1}^R \underbrace{\left(({\boldsymbol{u}}(t) \circ \boldsymbol{F}_{k} \circ \boldsymbol{\hat{\gamma}}_{n}(\hat{q}_r)) \cdot \boldsymbol{\nu}_{kn} \right)}_{\eqqcolon [{\boldsymbol{U}}_{\boldsymbol{\nu}}]_{k,n,r}} \, \underbrace{\omega_{r} \,  \hat{\varphi}_{i} \circ \boldsymbol{\hat{\gamma}}_{n}(\hat{q}_r) \, \hat{\mu}_{j} \circ \hat{\beta}_{kn}(\hat{q}_r)}_{\eqqcolon [\mathsf{\hat{S}}]_{i,j,n,r,l}} \,,
%> @f}
%> 
%> @f[
%> \mathsf{S} = \sum_{n=1}^3 \underbrace{\begin{bmatrix} 
%>   \delta_{E_{1n}=E_{1}} & \dots & \delta_{E_{1n}=E_{\bar{K}}} \\
%>   \vdots & \ddots & \vdots \\
%>   \delta_{E_{Kn}=E_{1}} & \dots & \delta_{E_{Kn}=E_{\bar{K}}}
%> \end{bmatrix}}_{\eqqcolon \mathsf{\Delta}_{n}} \otimes_\mathrm{V} \left( 
%> \sum_{r=1}^R \sum_{l=1}^2 \begin{bmatrix}
%>   \delta_{E_{1n}\in\ensuremath{\mathcal{E}}_{\text{int}}\cup\ensuremath{\mathcal{E}}_{\text{out}}} \, \ensuremath{|E_{1n}|} \, 
%>     \delta_{\kappa(\rho(1,n),l) = 1} \, [{\boldsymbol{U}}_{\boldsymbol{\nu}}]_{1,n,r} \\
%>   \vdots \\
%>   \delta_{E_{Kn}\in\ensuremath{\mathcal{E}}_{\text{int}}\cup\ensuremath{\mathcal{E}}_{\text{out}}} \, \ensuremath{|E_{Kn}|} \, 
%>     \delta_{\kappa(\rho(K,n),l) = K} \, [{\boldsymbol{U}}_{\boldsymbol{\nu}}]_{1,n,r} 
%> \end{bmatrix} \otimes [\mathsf{\hat{S}}]_{:,:,n,r,l} \right)\,.
%> @f]
%> 
%> 
%> All other entries are zero.
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param refEdgePhiIntMuPerQuad TODO
%> @param valOnQuad TODO
%>                    @f$[\bar{N} \times \bar{N}]@f$
%> @retval ret        The assembled matrix @f$[\bar{K}\bar{N} \times \bar{K}\bar{N}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
function ret = assembleMatEdgePhiIntMuVal(g, markE0T, refEdgePhiIntMuPerQuad, valOnQuad)
% Extract dimensions
K = g.numT;
bar{k} = g.numE;
[N, Nmu, ~, R] = size(refEdgePhiIntMuPerQuad{1});

% Check function arguments
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMuPerQuad, {'cell'}, {'size', [2 1]}, mfilename, 'refEdgePhiIntMuPerQuad');
validateattributes(refEdgePhiIntMuPerQuad{1}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'refEdgePhiIntMuPerQuad{1}');
validateattributes(refEdgePhiIntMuPerQuad{2}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'refEdgePhiIntMuPerQuad{2}');
validateattributes(valOnQuad, {'numeric'}, {'size', [K 3 R]}, mfilename, 'valOnQuad');

% Assemble matrix
ret = sparse(K*N, bar{k}*Nmu);
for n = 1 : 3
  RknTimesVal = sparse(K*N, Nmu);
  markAreaE0T = markE0T(:, n) .* g.areaE0T(:, n);
  for l = 1 : 2
    markAreaSideE0T = markAreaE0T .*  g.markSideE0T(:, n, l);
    for r = 1 : R
      RknTimesVal = RknTimesVal + kron(markAreaSideE0T .* valOnQuad(:, n, r), refEdgePhiIntMuPerQuad{l}(:, :, n, r));
    end % for r
  end % for l
  ret = ret + kronVec(sparse(1:K, g.E0T(:, n), ones(K, 1), K, bar{k}), RknTimesVal);
end % for n
end % function
