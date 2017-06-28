% Assembles two matrices, each containing integrals of products of a basis 
% function with a (spatial) derivative of a basis function and with a 
% discontinuous coefficient function.

%===============================================================================
%> @file assembleMatElemTetraDphiPhiFuncDisc.m
%>
%> @brief Assembles two matrices, each containing integrals of products of a 
%>        basis function with a (spatial) derivative of a basis function
%>        and with a discontinuous coefficient function.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{G}^m, m \in \{1,2\}@f$
%>        containing integrals of products of a basis function with a (spatial)
%>        derivative of a basis function and with a discontinuous coefficient 
%>        function.
%>
%> The matrices @f$\mathsf{G}^m \in \mathbb{R}^{KN\times KN}@f$ are block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{G}^m]_{(k-1)N+i,(k-1)N+j} = \sum_{l=1}^N D_{kl}(t) \int_{T_k} 
%>     \partial_{x^m}\varphi_{ki}\varphi_{kl}\varphi_{kj}\mathrm{d}\mathbf{x}\,.
%> @f]
%> All other entries are zero.
%> For the implementation, the element integrals are backtransformed to the
%> reference square @f$\hat{T} = [0,1]^2@f$ using a mapping
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhiPhi Local matrix @f$\hat{\mathsf{G}}@f$ as provided
%>                    by <code>integrateRefElemTetraDphiPhiPhi()</code>.
%>                    @f$[N \times N \times N \times 2]@f$
%> @param dataDisc    A representation of the discrete function ,e.g., as 
%>                    computed by <code>projectFuncCont2DataDiscTetra()</code>
%>                    @f$[K \times N]@f$
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
function ret = assembleMatElemTetraDphiPhiFuncDisc(g, refElemDphiPhiPhi, dataDisc)
if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    ret = assembleMatElemTetraDphiPhiFuncDiscMatrix(g, refElemDphiPhiPhi, dataDisc);
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    ret = assembleMatElemTetraDphiPhiFuncDiscVector(g, refElemDphiPhiPhi, dataDisc);
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  ret = assembleMatElemTetraDphiPhiFuncDiscScalar(g, refElemDphiPhiPhi, dataDisc);
end % if
end % function

function ret = assembleMatElemTetraDphiPhiFuncDiscMatrix(g, refElemDphiPhiPhi, dataDisc)
[K, N] = size(dataDisc{1,1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for r = 1 : 2
    for l = 1 : N
      for s = 1 : 3
        ret{m} = ret{m} + ...
                   kron(spdiags(dataDisc{r,m}(:,l) .* g.J0T{s}(:,3-r,3-r), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,  r)) - ...
                   kron(spdiags(dataDisc{r,m}(:,l) .* g.J0T{s}(:,3-r,  r), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,3-r));
      end % for s
    end % for l
  end % for r
end % for m
end  % function

function ret = assembleMatElemTetraDphiPhiFuncDiscVector(g, refElemDphiPhiPhi, dataDisc)
[K, N] = size(dataDisc{1});
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for l = 1 : N
    for s = 1 : 3
      ret{m} = ret{m} + ...
                 kron(spdiags(dataDisc{m}(:,l) .* g.J0T{s}(:,3-m,3-m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,  m)) - ...
                 kron(spdiags(dataDisc{m}(:,l) .* g.J0T{s}(:,3-m,  m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,3-m));
    end % for s
  end % for l
end % for m
end  % function

function ret = assembleMatElemTetraDphiPhiFuncDiscScalar(g, refElemDphiPhiPhi, dataDisc)
[K, N] = size(dataDisc);
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  for l = 1 : N
    for s = 1 : 3
      ret{m} = ret{m} + ...
                 kron(spdiags(dataDisc(:,l) .* g.J0T{s}(:,3-m,3-m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,  m)) - ...
                 kron(spdiags(dataDisc(:,l) .* g.J0T{s}(:,3-m,  m), 0, K, K), refElemDphiPhiPhi{s}(:,:,l,3-m));
    end % for s
  end % for l
end % for m
end  % function