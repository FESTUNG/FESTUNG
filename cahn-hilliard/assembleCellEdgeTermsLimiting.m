% Assembles a matrix containing integrals over interior edges of products of a 
% basis function, a derivative of a basis function, and the normal of the triangle.

%===============================================================================
%> @file 
%>
%> @brief Assembles a matrix containing integrals over interior edges of
%>        products of basis function, a derivative of a basis function,
%>        and the normal of the triangle.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals over interior edges of
%>        products of a basis function, a derivative of a basis function,
%>        and the normal of the triangle.
%>
%> The matrices @f$\mathsf{{U}} \in \mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal 
%> blocks. Diagonal blocks are defined as 
%> @f[
%> [\mathsf{{T}}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \int_{E_{kn}} \varphi_{ki} \nabla \varphi_{kj} \cdot \nu_{e_{kn}} \mathrm{d}s \,,
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{U}}]_{(k^--1)N+i,(k^+-1)N+j} =
%>  -\int_{E_{k^-n^-}} \varphi_{k^-i} \nabla \varphi_{k^+j} \cdot \nu_{e_{k^-n^-}} \mathrm{d}s \,.
%> @f]
%> Entries in off-diagonal blocks are potentially non-zero for pairs of elements
%> @f$T_{k^-}, T_{k^+}@f$ with @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$.
%> Note that the local edge index @f$n^-@f$ is given implicitly, since 
%> @f$\partial T_{k^-} \cap T_{k^+} \ne\emptyset@f$ consist of exactly one
%> edge @f$E_{k^-n^-} = E_{k^+n^+}@f$; all possible combinations of @f$n^-@f$ and @f$n^+@f$
%> must be considered.
%>
%> To allow for vectorization, the assembly is reformulated as described in
%> [FrankKuzminRupp2018].
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g., 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param refEdgeDPhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{U}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeDPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges}]@f$
%> @param refEdgeDPhiIntPhiExt Local matrix 
%>                    @f$\hat{\mathsf{{U}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefDEdgePhiIntPhiExt()</code>.
%>                    @f$[N \times N \times n_\mathrm{edges} \times n_\mathrm{edges}]@f$
%>                    @f$[N \times N \times n_\mathrm{edges} \times n_\mathrm{edges}]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Andreas Rupp, 2018.
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
function ret = assembleCellEdgeTermsLimiting(g, basesOnQuad, dataDisc, func, refEdgeDPhiIntPhiIntQuad, ...
    refEdgeDPhiIntPhiExtQuad, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, eta, sigma)

[K, N] = size(dataDisc); nEdges = size(g.E0T, 2);
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);
[~, W] = quadRule1D(qOrd); R = length(W);

funcValEdge = zeros(K,size(basesOnQuad.phi1D{qOrd}(:, 1),1),3);
for nn = 1 : 3
  for l = 1 : N 
    funcValEdge(:,:,nn) = funcValEdge(:,:,nn) + kron(dataDisc(:,l) , basesOnQuad.phi1D{qOrd}(:, l, nn).');
  end %for l = 1 : N
end %for nn = 1 : 3
funcValEdge = func(funcValEdge);

% Create the transposed of Tdiag
refEdgePhiIntDPhiIntQuad = cell(2,1);
refEdgePhiIntDPhiIntQuad{1} = permute(refEdgeDPhiIntPhiIntQuad{1}, [2 1 3 4]);
refEdgePhiIntDPhiIntQuad{2} = permute(refEdgeDPhiIntPhiIntQuad{2}, [2 1 3 4]);

% Create the transposed of Toffdiag (where also the n-indices are switched)
refEdgePhiIntDPhiExtQuad = cell(2,1);
refEdgePhiIntDPhiExtQuad{1} = permute(refEdgeDPhiIntPhiExtQuad{1}, [2 1 4 3 5]);
refEdgePhiIntDPhiExtQuad{2} = permute(refEdgeDPhiIntPhiExtQuad{2}, [2 1 4 3 5]);

markE0T = g.markE0Tint;

ret = cell(nEdges,1);
for i = 1 : nEdges
  ret{i} = sparse(K, K*N);
end
indexVecHelp = (1 : K).';

if numel(g.markE0TE0T) == nEdges % mapping from nn to np implicitly given

  assert(0 == 1, 'Not yet implemented!');
  
else % mapping from nn to np explicitly given
  
  for r = 1 : R
    for nn = 1 : nEdges
      markCoefE0T = markE0T(:,nn) .* g.areaE0T(:,nn) ./ (2 * g.areaT);
      ret{nn} = ret{nn} ...
        - 0.5 * kron(spdiags(funcValEdge(:,r,nn) .* (markCoefE0T .* ( g.B(:,2,2) .* g.nuE0T(:,nn,1) - g.B(:,1,2) .* g.nuE0T(:,nn,2) )) , 0, K, K), refEdgePhiIntDPhiIntQuad{1}(1,:,nn,r)) ...
        - 0.5 * kron(spdiags(funcValEdge(:,r,nn) .* (markCoefE0T .* (-g.B(:,2,1) .* g.nuE0T(:,nn,1) + g.B(:,1,1) .* g.nuE0T(:,nn,2) )) , 0, K, K), refEdgePhiIntDPhiIntQuad{2}(1,:,nn,r));
      for np = 1 : nEdges
        if nnz(g.markE0TE0T{nn, np}) > 0
          indexVec = g.markE0TE0T{nn, np} * indexVecHelp;
          indexVec(indexVec == 0) = 1;
          markCoefE0Tind = markE0T(:,nn) .* g.areaE0T(:,nn) ./ (2 * g.areaT(indexVec));
          ret{nn} = ret{nn} ...
            - 0.5 * kron(bsxfun(@times, g.markE0TE0T{nn, np}, funcValEdge(indexVec,r, np) .* markCoefE0Tind .* ( g.B(indexVec,2,2) .* g.nuE0T(:,nn,1) - g.B(indexVec,1,2) .* g.nuE0T(:,nn,2) )), refEdgePhiIntDPhiExtQuad{1}(1,:,nn,np,r)) ...
            - 0.5 * kron(bsxfun(@times, g.markE0TE0T{nn, np}, funcValEdge(indexVec,r, np) .* markCoefE0Tind .* (-g.B(indexVec,2,1) .* g.nuE0T(:,nn,1) + g.B(indexVec,1,1) .* g.nuE0T(:,nn,2) )), refEdgePhiIntDPhiExtQuad{2}(1,:,nn,np,r));
        end % for nnz > 0
      end % for np
    end % for nn
  end % for r
  
  for r = 1 : R
    for nn = 1 : nEdges
      markCoefE0T = markE0T(:,nn) .* g.areaE0T(:,nn) ./ (2 * g.areaT);
      ret{nn} = ret{nn} ...
        + 0.5 * eta * kron(spdiags(funcValEdge(:,r,nn) .* markCoefE0T .*  (g.B(:,2,2) .* g.nuE0T(:,nn,1) - g.B(:,1,2) .* g.nuE0T(:,nn,2) ) , 0, K, K), refEdgeDPhiIntPhiIntQuad{1}(1,:,nn,r)) ...
        + 0.5 * eta * kron(spdiags(funcValEdge(:,r,nn) .* markCoefE0T .* (-g.B(:,2,1) .* g.nuE0T(:,nn,1) + g.B(:,1,1) .* g.nuE0T(:,nn,2) ) , 0, K, K), refEdgeDPhiIntPhiIntQuad{2}(1,:,nn,r));
      for np = 1 : nEdges
        if nnz(g.markE0TE0T{nn, np}) > 0
          ret{nn} = ret{nn} ...
            - 0.5 * eta * kron(bsxfun(@times, g.markE0TE0T{nn, np}, funcValEdge(:,r,nn) .* markCoefE0T .* ( g.B(:,2,2) .* g.nuE0T(:,nn,1) - g.B(:,1,2) .* g.nuE0T(:,nn,2) )), refEdgeDPhiIntPhiExtQuad{1}(1,:,nn,np,r)) ...
            - 0.5 * eta * kron(bsxfun(@times, g.markE0TE0T{nn, np}, funcValEdge(:,r,nn) .* markCoefE0T .* (-g.B(:,2,1) .* g.nuE0T(:,nn,1) + g.B(:,1,1) .* g.nuE0T(:,nn,2) )), refEdgeDPhiIntPhiExtQuad{2}(1,:,nn,np,r));
        end
      end % for np
    end % for nn
  end % for r
  
   for nn = 1 : nEdges
    markCoefE0T = markE0T(:,nn);
    ret{nn} = ret{nn} + sigma * kron(spdiags(markCoefE0T, 0, K, K), refEdgePhiIntPhiInt(1,:,nn));
    for np = 1 : nEdges
      ret{nn} = ret{nn} - sigma * kron(bsxfun(@times, g.markE0TE0T{nn, np}, markCoefE0T), refEdgePhiIntPhiExt(1,:,nn,np));
    end % for np
  end % for nn
  
end % if mapping for nn
end % function
