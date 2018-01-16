% Assembles two matrices containing integrals over edges of products of two 
% basis functions from the interior of each element and a function in discrete
% representation with a component of the edge normal.

%===============================================================================
%> @file assembleMatEdgePhiIntPhiIntFuncDiscIntNu.m
%>
%> @brief Assembles two matrices containing integrals over edges of products of
%>        two basis functions from the interior of each element and a function
%>        in discrete representation with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles matrices
%>        @f$\mathsf{{R}}^m_\mathrm{D}, m\in\{1,2\}@f$ containing integrals over  
%>        edges of products of two basis functions and a function in discrete 
%>        representation from the interior of each element.
%>
%> The matrix @f$\mathsf{{R}}^m_\mathrm{D}\in\mathbb{R}^{KN\times KN}@f$
%> is block diagonal and defined as
%> @f[
%> [\mathsf{{R}}^m_\mathrm{D}]_{(k-1)N+i,(k-1)N+j} =
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_D}
%>  \nu_{kn}^m \sum_{l=1}^{N_D} D_{kl}(t) 
%>  \int_{E_kn} \varphi_{ki}\varphi_{kl}\varphi_{kj} \,,
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the unit
%> normal and @f$D_{kl}(t)@f$ the coefficients of the discrete representation
%> of a function 
%> @f$ d_h(t, \mathbf{x}) = \sum_{l=1}^{N_D} D_{kl}(t) \varphi_{kl}(\mathbf{x}) @f$
%> on an element @f$T_k@f$, where @f$N_D@f$ is the number of local degrees of freedom
%> of the space in which @f$d_h(t,.)@f$ is projected.
%> All other entries are zero.
%>
%> To allow for vectorization, the assembly is reformulated as
%> @f[
%> \mathsf{{R}}^m_\mathrm{D} = \sum_{n=1}^{n_\mathrm{edges}} \sum_{l=1}^{N_D}
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\mathrm{D}} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\mathrm{D}}
%>   \end{bmatrix} \circ
%>   \begin{bmatrix}
%>     \nu^m_{1n} |E_{1n}| D_{1l}(t) & & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} |E_{Kn}| D_{Kl}(t)
%>   \end{bmatrix} \otimes [\hat{\mathsf{{R}}}^\mathrm{diag}]_{:,:,l,n}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\mathrm{D}}@f$ denotes the Kronecker 
%> delta, @f$\circ@f$ denotes the Hadamard product, and @f$\otimes@f$ denotes
%> the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{R}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times N_D\times{n_\mathrm{edges}}}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{diag}]_{i,j,l,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_l \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>.
%>
%> It is essentially the same as the diagonal part of
%> <code>assembleMatEdgePhiPhiFuncDiscNu()</code>.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times {n_\mathrm{edges}}]@f$
%> @param refEdgePhiIntPhiIntPhiInt  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiIntPhiInt()</code>.
%>                    @f$[N \times N \times N \times {n_\mathrm{edges}}]@f$
%> @param dataDisc    A representation of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$
%> @param areaNuMarkE0T (optional) argument to provide precomputed values
%>                    for the products of <code>markE0T</code>,
%>                    <code>g.areaE0T</code>, and <code>g.nuE0T</code>
%>                    @f$[2 \times 1 \text{ cell}]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
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
function ret = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, dataDisc, areaNuMarkE0T)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

% Select function depending on type of dataDisc
if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    dataN = size(dataDisc{1, 1}, 2);
    fn = @assembleMatEdgePhiIntPhiIntFuncDiscIntMatrixNu;
    
    validateattributes(dataDisc{1, 1}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{1, 1}');
    validateattributes(dataDisc{1, 2}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{1, 2}');
    validateattributes(dataDisc{2, 1}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{2, 1}');
    validateattributes(dataDisc{2, 2}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{2, 2}');
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    dataN = size(dataDisc{1, 1}, 2);
    fn = @assembleMatEdgePhiIntPhiIntFuncDiscIntVectorNu;
    
    validateattributes(dataDisc{1}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{1}');
    validateattributes(dataDisc{2}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{2}');
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  dataN = size(dataDisc, 2);
  fn = @assembleMatEdgePhiIntPhiIntFuncDiscIntScalarNu;
  
  validateattributes(dataDisc, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc');
end % if

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntPhiIntPhiInt, {'numeric'}, {'size', [N N dataN nEdges]}, mfilename, 'refEdgePhiIntPhiIntPhiInt');

if nargin < 5
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

ret = fn(g, refEdgePhiIntPhiIntPhiInt, dataDisc, areaNuMarkE0T);
end % function

function ret = assembleMatEdgePhiIntPhiIntFuncDiscIntMatrixNu(g, refEdgePhiIntPhiIntPhiInt, dataDisc, areaNuMarkE0T)
[K, dataN] = size(dataDisc{1, 1});  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  Rm = zeros(K*N, N);
  for n = 1 : nEdges
    for r = 1 : 2
      for l = 1 : dataN
        Rm = Rm + kron(areaNuMarkE0T{r}(:, n) .* dataDisc{r, m}(:, l), refEdgePhiIntPhiIntPhiInt(:, :, l, n));
      end % for l
    end % for r
  end % for n
  ret{m} = kronVec(speye(K, K), Rm);
end  % for m
end % function

function ret = assembleMatEdgePhiIntPhiIntFuncDiscIntVectorNu(g, refEdgePhiIntPhiIntPhiInt, dataDisc, areaNuMarkE0T)
[K, dataN] = size(dataDisc{1});  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  Rm = zeros(K*N, N);
  for n = 1 : nEdges
    for l = 1 : dataN
      Rm = Rm + kron(areaNuMarkE0T{m}(:, n) .* dataDisc{m}(:, l), refEdgePhiIntPhiIntPhiInt(:, :, l, n));
    end % for l
  end % for n
  ret{m} = kronVec(speye(K, K), Rm);
end  % for m
end % function

function ret = assembleMatEdgePhiIntPhiIntFuncDiscIntScalarNu(g, refEdgePhiIntPhiIntPhiInt, dataDisc, areaNuMarkE0T)
[K, dataN] = size(dataDisc);  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for m = 1 : 2
  Rm = zeros(K*N, N);
  for n = 1 : nEdges
    for l = 1 : dataN
      Rm = Rm + kron(areaNuMarkE0T{m}(:, n) .* dataDisc(:, l), refEdgePhiIntPhiIntPhiInt(:, :, l, n));
    end % for l
  end % for n
  ret{m} = kronVec(speye(K, K), Rm);
end  % for m
end % function
