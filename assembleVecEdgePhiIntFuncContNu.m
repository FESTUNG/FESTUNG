% Assembles two vectors containing integrals over edges of products of a basis 
% function with a continuous function and one of the components of the edge 
% normal.

%===============================================================================
%> @file
%>
%> @brief Assembles two vectors containing integrals over edges of products of a 
%>        basis function with a continuous function and one of the components of
%>        the edge normal.
%===============================================================================
%>
%> @brief Assembles two vectors @f$\mathbf{J}_\mathrm{D}^m, m\in\{1,2\}@f$ 
%>        containing integrals over edges of products of a basis function with a
%>        continuous function @f$c_\mathrm{D}(t, \mathbf{x})@f$ and one of the 
%>        components of the edge normal.
%>
%> The vectors @f$\mathbf{J}_\mathrm{D}^m \in \mathbb{R}^{KN}@f$ are defined
%> component-wise by
%> @f[
%> [\mathbf{J}_\mathrm{D}^m]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \nu_{kn}^m \int_{E_{kn}} \varphi_{ki} c_\mathrm{D}(t) \mathrm{d}s\,,
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the edge
%> normal.
%>
%> For the implementation, the integrals are backtransformed to the
%> reference element @f$\hat{T}@f$ and further to the unit interval @f$[0,1]@f$
%> using the mapping @f$\hat{\mathbf{\gamma}}_n(s)@f$ as provided by 
%> <code>gammaMap()</code> or <code>gammaMapTetra()</code>.
%> This gives the component-wise formulation
%> @f[
%>  [\mathbf{J}_\mathrm{D}^m]_{(k-1)N+i} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} \nu_{kn}^m |E_{kn}|
%>  \int_0^1 \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(s)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(s))
%>     \mathrm{d}s \,.
%> @f]
%> This integral is then approximated using a 1D quadrature rule provided by
%> <code>quadRule1D()</code> 
%> @f[
%>  [\mathbf{J}_\mathrm{D}^m]_{(k-1)N+i} \approx
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D} \nu_{kn}^m |E_{kn}|
%>  \sum_{r=1}^R \omega_r \hat{\varphi}_{i} \circ \hat{\mathbf{\gamma}}_n(q_r)
%>     c_\mathrm{D}(t, \mathbf{F}_k \circ \hat{\mathbf{\gamma}}_n(q_r)) \,,
%> @f]
%> allowing to vectorize over all triangles.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the vector entries should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param  funcCont   A function handle for the continuous function
%> @param  N          The number of local degrees of freedom @f$[\text{scalar}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                    functions on quadrature points. Must provide at
%>                    least phi1D.
%> @param  qOrd       (optional) The order of the quadrature rule to be used. 
%>                    Defaults to @f$2p+1 @f$.
%> @param  areaNuMarkE0T  (optional) argument to provide precomputed values
%>                    for the products of <code>markE0T</code>,
%>                    <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                    @f$[2 \times 1 \text{ cell}]@f$
%> @retval ret        The assembled vectors @f$2\times1 \text{ cell}@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2016.
%> @author Balthasar Reuter, 2017.
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
function ret = assembleVecEdgePhiIntFuncContNu(g, markE0T, funcCont, N, basesOnQuad, qOrd, areaNuMarkE0T)
% Extract dimensions
K = g.numT;  nEdges = size(g.E0T, 2);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

% Determine quadrature rule and mapping to physical element
if nargin < 6
  switch nEdges
    case 3
      p = (sqrt(8*N+1)-3)/2;
    case 4
      p = sqrt(N)-1;
    otherwise
      error('Unknown number of edges')
  end % switch
  qOrd = 2*p+1;  
end % if
[Q, W] = quadRule1D(qOrd);

% Determine helper array
if nargin < 7
  areaNuMarkE0T = cell(2,1);
  if isfield(g, 'areaNuE0T')
    areaNuMarkE0T{1} = markE0T .* g.areaNuE0T(:, :, 1);
    areaNuMarkE0T{2} = markE0T .* g.areaNuE0T(:, :, 2);
  else
    areaNuMarkE0T{1} = markE0T .* g.areaE0T .* g.nuE0T(:, :, 1);
    areaNuMarkE0T{2} = markE0T .* g.areaE0T .* g.nuE0T(:, :, 2);
  end % if
end % if

validateattributes(areaNuMarkE0T, {'cell'}, {'numel', 2}, mfilename, 'areaNuMarkE0T');
validateattributes(areaNuMarkE0T{1}, {'numeric'}, {'size', [K nEdges]}, mfilename, 'areaNuMarkE0T{1}');
validateattributes(areaNuMarkE0T{2}, {'numeric'}, {'size', [K nEdges]}, mfilename, 'areaNuMarkE0T{2}');

% Determine edge mapping
switch nEdges
  case 3
    gamma = @gammaMap;
  case 4
    gamma = @gammaMapTetra;
  otherwise
    error('Unknown number of edges')
end % switch

% Assemble vector
ret = { zeros(K, N); zeros(K, N) };

if iscell(funcCont)
  validateattributes(funcCont, {'cell'}, {'numel', 2}, mfilename, 'funcCont');
  validateattributes(funcCont{1}, {'function_handle'}, {}, mfilename, 'funcCont{1}');
  validateattributes(funcCont{2}, {'function_handle'}, {}, mfilename, 'funcCont{2}');
  
  for n = 1 : nEdges
    [Q1, Q2] = gamma(n, Q);
    X1 = g.mapRef2Phy(1, Q1, Q2);
    X2 = g.mapRef2Phy(2, Q1, Q2);
    for m = 1 : 2
      funcQ0E = funcCont{m}(X1, X2);
      for i = 1 : N
        integral = funcQ0E * ( W' .* basesOnQuad.phi1D{qOrd}(:, i, n) );
        ret{m}(:, i) = ret{m}(:, i) + areaNuMarkE0T{m}(:, n) .* integral;
      end % for i
    end % for m
  end % for n
else
  validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
  
  for n = 1 : nEdges
    [Q1, Q2] = gamma(n, Q);
    funcQ0E = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
    for i = 1 : N
      integral = funcQ0E * ( W' .* basesOnQuad.phi1D{qOrd}(:, i, n) );
      ret{1}(:, i) = ret{1}(:, i) + areaNuMarkE0T{1}(:, n) .* integral;
      ret{2}(:, i) = ret{2}(:, i) + areaNuMarkE0T{2}(:, n) .* integral;
    end % for i
  end % for n
end % if

ret{1} = reshape(ret{1}', K*N, 1);  
ret{2} = reshape(ret{2}', K*N, 1);
end % function
% 
% 
% if nargin > 5
%   ret = assembleVecEdgePhiIntFuncContNu_withAreaNuE0Tbdr(g, funcCont, N, basesOnQuad, markAreaNuE0T);
% elseif isfield(g, 'areaNuE0T')
%   ret = assembleVecEdgePhiIntFuncContNu_noAreaNuE0Tbdr_withAreaNuE0T(g, markE0T, funcCont, N, basesOnQuad);
% else
%   ret = assembleVecEdgePhiIntFuncContNu_noAreaNuE0Tbdr_noAreaNuE0T(g, markE0T, funcCont, N, basesOnQuad);
% end % if
% 
% end % function
% %
% %===============================================================================
% %> @brief Helper function for the case that assembleVecEdgePhiIntFuncContNu()
% %> was called with a precomputed field areaNuE0Tbdr.
% %
% function ret = assembleVecEdgePhiIntFuncContNu_withAreaNuE0Tbdr(g, funcCont, N, basesOnQuad, areaNuE0Tbdr)
% % Determine quadrature rule and mapping to physical element
% K = g.numT; p = (sqrt(8*N+1)-3)/2;
% qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd); 
% 
% % Assemble vector
% ret = cell(2, 1);  
% ret{1} = zeros(K, N);  
% ret{2} = zeros(K, N);
% for n = 1 : 3
%   [Q1, Q2] = gammaMap(n, Q);
%   cDn = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
%   for i = 1 : N
%     integral = cDn * ( W .* basesOnQuad.phi1D{qOrd}(:, i, n)' )';
%     ret{1}(:,i) = ret{1}(:,i) + areaNuE0Tbdr{n,1} .* integral;
%     ret{2}(:,i) = ret{2}(:,i) + areaNuE0Tbdr{n,2} .* integral;
%   end % for
% end % for
% 
% ret{1} = reshape(ret{1}',K*N,1);  
% ret{2} = reshape(ret{2}',K*N,1);
% end % function
% %
% %===============================================================================
% %> @brief Helper function for the case that assembleVecEdgePhiIntFuncContNu()
% %> was called without a precomputed field areaNuE0Tbdr but parameter g provides
% %> a field areaNuE0T.
% %
% function ret = assembleVecEdgePhiIntFuncContNu_noAreaNuE0Tbdr_withAreaNuE0T(g, markE0Tbdr, funcCont, N, basesOnQuad)
% % Determine quadrature rule and mapping to physical element
% K = g.numT; p = (sqrt(8*N+1)-3)/2;
% qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd); 
% 
% % Assemble vector
% ret = cell(2, 1);  
% ret{1} = zeros(K, N);  
% ret{2} = zeros(K, N);
% for n = 1 : 3
%   [Q1, Q2] = gammaMap(n, Q);
%   cDn = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
%   for i = 1 : N
%     integral = cDn * ( W .* basesOnQuad.phi1D{qOrd}(:, i, n)' )';
%     ret{1}(:,i) = ret{1}(:,i) + markE0Tbdr(:,n).*g.areaNuE0T{n,1} .* integral;
%     ret{2}(:,i) = ret{2}(:,i) + markE0Tbdr(:,n).*g.areaNuE0T{n,2} .* integral;
%   end % for
% end % for
% 
% ret{1} = reshape(ret{1}',K*N,1);  
% ret{2} = reshape(ret{2}',K*N,1);
% end % function
% %
% %===============================================================================
% %> @brief Helper function for the case that assembleVecEdgePhiIntFuncContNu()
% %> was called without a precomputed field areaNuE0Tbdr and parameter g provides
% %> no field areaNuE0T.
% %
% function ret = assembleVecEdgePhiIntFuncContNu_noAreaNuE0Tbdr_noAreaNuE0T(g, markE0Tbdr, funcCont, N, basesOnQuad)
% % Determine quadrature rule and mapping to physical element
% K = g.numT; p = (sqrt(8*N+1)-3)/2;
% qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd); 
% 
% % Assemble vector
% ret = cell(2, 1);  
% ret{1} = zeros(K, N); 
% ret{2} = zeros(K, N);
% for n = 1 : 3
%   [Q1, Q2] = gammaMap(n, Q);
%   cDn = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
%   Jkn = markE0Tbdr(:,n) .* g.areaE0T(:,n);
%   for i = 1 : N
%     integral = cDn * ( W .* basesOnQuad.phi1D{qOrd}(:, i, n)' )';
%     ret{1}(:,i) = ret{1}(:,i) + Jkn .* g.nuE0T(:,n,1) .* integral;
%     ret{2}(:,i) = ret{2}(:,i) + Jkn .* g.nuE0T(:,n,2) .* integral;
%   end % for
% end % for
% 
% ret{1} = reshape(ret{1}',K*N,1);  
% ret{2} = reshape(ret{2}',K*N,1);
% end % function