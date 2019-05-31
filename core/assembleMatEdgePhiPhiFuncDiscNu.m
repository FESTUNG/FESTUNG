% Assembles two matrices containing integrals over interior edges of products of
% two basis functions with a discontinuous coefficient function and with a
% component of the edge normal.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices containing integrals over interior edges of 
%>        products of two basis functions with a discontinuous coefficient 
%>        function and with a component of the edge normal.
%===============================================================================
%>
%> @brief Assembles the matrices @f$\mathsf{{R}}^m, m\in\{1,2\}@f$ containing 
%>        integrals over interior edges of products of two basis functions with 
%>        a discontinuous coefficient function and a component of the edge normal.
%>
%> The matrices @f$\mathsf{{R}}^m\in\mathbb{R}^{KN\times KN}@f$
%> consist of two kinds of contributions: diagonal blocks and off-diagonal
%> blocks. Diagonal blocks are defined as
%> @f[
%> [\mathsf{{R}}^m]_{(k-1)N+i,(k-1)N+j} = \frac{1}{2} 
%>  \sum_{E_{kn} \in \partial T_k \cap \mathcal{E}_\Omega}
%>  \nu_{kn}^m \sum_{l=1}^N D_{kl}(t) 
%>  \int_{E_kn} \varphi_{ki}\varphi_{kl}\varphi_{kj} \,,
%> @f]
%> and off-diagonal blocks are defined as
%> @f[
%> [\mathsf{{R}}^m]_{(k^--1)N+i,(k^+-1)N+j} =
%>   \frac{1}{2} \nu^m_{k^-n^-} \sum_{l=1}^N D_{k^+l}(t) \int_{E_{k^-n^-}} 
%>   \varphi_{k^-i} \varphi_{k^+l} \varphi_{k^+j} \mathrm{d}s \,.
%> @f]
%> with @f$\nu_{kn}@f$ the @f$m@f$-th component (@f$m\in\{1,2\}@f$) of the unit
%> normal and @f$D_{kl}(t)@f$ the coefficients of the discrete representation
%> of a function 
%> @f$ d_h(t, \mathbf{x}) = \sum_{l=1}^N D_{kl}(t) \varphi_{kl}(\mathbf{x}) @f$
%> on an element @f$T_k@f$.
%> All other entries are zero.
%>
%> To allow for vectorization, the assembly is reformulated as
%> @f$\mathsf{{R}}^m = \mathsf{{R}}^{m,\mathrm{diag}} + 
%>    \mathsf{{R}}^{m,\mathrm{offdiag}}@f$ with the blocks defined as
%> @f[
%> \mathsf{{R}}^{m,\mathrm{diag}} = \frac{1}{2} \sum_{n=1}^{n_\mathrm{edges}} \sum_{l=1}^N
%>   \begin{bmatrix}
%>     \delta_{E_{1n}\in\mathcal{E}_\Omega} &   & \\
%>     & ~\ddots~ & \\
%>     &          & \delta_{E_{Kn}\in\mathcal{E}_\Omega}
%>   \end{bmatrix} \circ
%>   \begin{bmatrix}
%>     \nu^m_{1n} |E_{1n}| D_{1l}(t) & & \\
%>     & ~\ddots~ & \\
%>     &          & \nu^m_{Kn} |E_{Kn}| D_{Kl}(t)
%>   \end{bmatrix} \otimes [\hat{\mathsf{{R}}}^\mathrm{diag}]_{:,:,l,n}\;,
%> @f]
%> and
%> @f[
%> \mathsf{{R}}^{m,\mathrm{offdiag}} = \frac{1}{2} 
%>   \sum_{n^-=1}^{n_\mathrm{edges}}\sum_{n^+=1}^{n_\mathrm{edges}} \sum_{l=1}^N
%>   \begin{bmatrix}
%>     0&\delta_{E_{1n^-} = E_{2n^+}}&\dots&\dots&\delta_{E_{1n^-}=E_{Kn^+}} \\
%>     \delta_{E_{2n^-} = E_{1n^+}}&0&\ddots& &\vdots \\
%>     \vdots & \ddots & \ddots & \ddots & \vdots \\
%>     \vdots & & \ddots & 0 & \delta_{E_{(K-1)n^-}=E_{Kn^+}} \\
%>     \delta_{E_{Kn^-} = E_{1n^+}}&\dots&\dots&\delta_{E_{Kn^-} = E_{(K-1)n^+}} &0
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     \nu_{1n^-}^m |E_{1n^-}| & \dots & \nu_{1n^-}^m |E_{1n^-}| \\
%>     \vdots & & \vdots \\
%>     \nu_{Kn^-}^m |E_{Kn^-}| & \dots & \nu_{Kn^-}^m |E_{Kn^-}|
%>   \end{bmatrix} \circ \begin{bmatrix}
%>     D_{1l}(t) & \dots & D_{Kl}(t) \\
%>     \vdots & & \vdots \\
%>     D_{1l}(t) & \dots & D_{Kl}(t)
%>   \end{bmatrix} \otimes 
%> [\hat{\mathsf{{R}}}^\mathrm{offdiag}]_{:,:,l,n^-,n^+}\;,
%> @f]
%> where @f$\delta_{E_{kn}\in\mathcal{E}_\Omega}, \delta_{E_{in^-}=E_{jn^+}}@f$ 
%> denote the Kronecker delta, @f$\circ@f$ denotes the Hadamard product, and 
%> @f$\otimes@f$ denotes the Kronecker product.
%>
%> The entries of matrix 
%> @f$\hat{\mathsf{{R}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times N\times{n_\mathrm{edges}}}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{diag}]_{i,j,l,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_l \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is defined in 
%> <code>gammaMap()</code>. The entries of matrix 
%> @f$\hat{\mathsf{{R}}}^\mathrm{offdiag} \in
%>    \mathbb{R}^{N\times N\times N\times {n_\mathrm{edges}}\times {n_\mathrm{edges}}}@f$
%> are given by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{offdiag}]_{i,j,l,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_l\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s)
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in
%> <code>theta()</code>.
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
%>                    @f$[N \times N \times N_\mathrm{data} \times {n_\mathrm{edges}}]@f$
%> @param refEdgePhiIntPhiExtPhiExt  Local matrix 
%>                    @f$\hat{\mathsf{{R}}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiExtPhiExt()</code>.
%>                    @f$[N \times N \times N_\mathrm{data} \times {n_\mathrm{edges}}]@f$ or
%>                    @f$[N \times N \times N_\mathrm{data} \times {n_\mathrm{edges}} \times {n_\mathrm{edges}}]@f$
%> @param dataDisc    A representation of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N_\mathrm{data}]@f$
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
function ret = assembleMatEdgePhiPhiFuncDiscNu(g, markE0T, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc, areaNuMarkE0T)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

% Select function depending on type of dataDisc
if iscell(dataDisc)
  if isequal(size(dataDisc), [2 2])
    dataN = size(dataDisc{1, 1}, 2);
    fn = @assembleMatEdgePhiPhiFuncDiscMatrixNu;
    
    validateattributes(dataDisc{1, 1}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{1, 1}');
    validateattributes(dataDisc{1, 2}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{1, 2}');
    validateattributes(dataDisc{2, 1}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{2, 1}');
    validateattributes(dataDisc{2, 2}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{2, 2}');
  elseif isvector(dataDisc) && numel(dataDisc) == 2
    dataN = size(dataDisc{1, 1}, 2);
    fn = @assembleMatEdgePhiPhiFuncDiscVectorNu;
    
    validateattributes(dataDisc{1}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{1}');
    validateattributes(dataDisc{2}, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc{2}');
  else
    error('dataDisc must be a KxN-matrix or a 2x1-cell or a 2x2-cell of such matrices.')
  end % if
else
  dataN = size(dataDisc, 2);
  fn = @assembleMatEdgePhiPhiFuncDiscScalarNu;
  
  validateattributes(dataDisc, {'numeric'}, {'size', [K dataN]}, mfilename, 'dataDisc');
end % if

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K nEdges]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntPhiIntPhiInt, {'numeric'}, {'size', [N N dataN nEdges]}, mfilename, 'refEdgePhiIntPhiIntPhiInt');

if nargin < 6
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

ret = fn(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc, areaNuMarkE0T);
end % function

function ret = assembleMatEdgePhiPhiFuncDiscMatrixNu(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc, areaNuMarkE0T)
[K, dataN] = size(dataDisc{1, 1});  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

ret = { sparse(K*N, K*N), sparse(K*N, K*N) };

if numel(g.markE0TE0T) == nEdges % mapping from nn to np implicitly given
  
  validateattributes(refEdgePhiIntPhiExtPhiExt, {'numeric'}, {'size', [N N N nEdges]}, mfilename, 'refEdgePhiIntPhiExtPhiExt');
    
  % Diagonal blocks
  for m = 1 : 2
    RtildeT = zeros(K*N, N);
    for n = 1 : nEdges
      for r = 1 : 2
        for l = 1 : dataN
          RtildeT = RtildeT + kron(0.5 * areaNuMarkE0T{r}(:, n) .* dataDisc{r, m}(:, l), refEdgePhiIntPhiIntPhiInt(:, :, l, n));
        end % for l
      end % for r
    end % for n
    ret{m} = kronVec(speye(K, K), RtildeT);
  end % for m
    
  % Off-diagonal blocks
  for n = 1 : nEdges
    for r = 1 : 2
      markOffdiag = (spdiags(0.5 * areaNuMarkE0T{r}(:, n), 0, K, K) * g.markE0TE0T{n}).';
      for m = 1 : 2
        % Coefficient-function for off-diagonal blocks
        RtildeT = zeros(K*N, N);
        for l = 1 : dataN
          RtildeT = RtildeT + kron(dataDisc{r, m}(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, n).');
        end % for l
        % Off-diagonal blocks
        ret{m} = ret{m} + kronVec(markOffdiag, RtildeT).';
      end % for m
    end % for r
  end  % for n
  
else % mapping from nn to np explicitly given
  
  validateattributes(refEdgePhiIntPhiExtPhiExt, {'numeric'}, {'size', [N N N nEdges nEdges]}, mfilename, 'refEdgePhiIntPhiExtPhiExt');

  for nn = 1 : nEdges
    for r = 1 : 2
      for m = 1 : 2
        % Off-diagonal blocks
        for np = 1 : nEdges
          RtildeT = zeros(K*N, N);
          for l = 1 : dataN
            RtildeT = RtildeT + kron(dataDisc{r, m}(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, nn, np).');
          end % for l
          if isfield(g, 'areaNuE0TE0T')
            ret{m} = ret{m} + kronVec(0.5 * g.areaNuE0TE0T{nn, np, r}.', RtildeT).';
          else
            ret{m} = ret{m} + kronVec(0.5 * bsxfun(@times, g.markE0TE0T{nn, np}, areaNuMarkE0T{r}(:, nn)).', RtildeT).';
          end % if
        end % for np
        % Diagonal blocks
        for l = 1 : dataN
          ret{m} = ret{m} + kron(spdiags(0.5 * areaNuMarkE0T{r}(:, nn) .* dataDisc{r, m}(:, l), 0, K, K), refEdgePhiIntPhiIntPhiInt(:, :, l, nn));
        end % for l
      end % for m
    end % for r
  end % for nn
  
end % if
end % function

function ret = assembleMatEdgePhiPhiFuncDiscVectorNu(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc, areaNuMarkE0T)
[K, dataN] = size(dataDisc{1});  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

ret = { sparse(K*N, K*N), sparse(K*N, K*N) };

if numel(g.markE0TE0T) == nEdges % mapping from nn to np implicitly given
  
  validateattributes(refEdgePhiIntPhiExtPhiExt, {'numeric'}, {'size', [N N N nEdges]}, mfilename, 'refEdgePhiIntPhiExtPhiExt');
  
  % Diagonal blocks
  for m = 1 : 2
    RtildeT = zeros(K*N, N);
    for n = 1 : nEdges
      for l = 1 : dataN
        RtildeT = RtildeT + kron(0.5 * areaNuMarkE0T{m}(:, n) .* dataDisc{m}(:, l), refEdgePhiIntPhiIntPhiInt(:, :, l, n));
      end % for l
    end % for n
    ret{m} = kronVec(speye(K, K), RtildeT);
  end % for m
    
  % Off-diagonal blocks
  for n = 1 : nEdges
    for m = 1 : 2
      % Coefficient-function for off-diagonal blocks
      RtildeT = zeros(K*N, N);
      for l = 1 : dataN
        RtildeT = RtildeT + kron(dataDisc{m}(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, n).');
      end % for l
      % Off-diagonal blocks
      ret{m} = ret{m} + kronVec((spdiags(0.5 * areaNuMarkE0T{m}(:, n), 0, K, K) * g.markE0TE0T{n}).', RtildeT).';
    end % for m
  end  % for n
  
  
else % mapping from nn to np explicitly given
  
  validateattributes(refEdgePhiIntPhiExtPhiExt, {'numeric'}, {'size', [N N N nEdges nEdges]}, mfilename, 'refEdgePhiIntPhiExtPhiExt');

  for nn = 1 : nEdges
    for m = 1 : 2
      % Off-diagonal blocks
      for np = 1 : nEdges
        RtildeT = zeros(K*N, N);
        for l = 1 : dataN
          RtildeT = RtildeT + kron(dataDisc{m}(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, nn, np).');
        end % for l
        if isfield(g, 'areaNuE0TE0T')
          ret{m} = ret{m} + kronVec(0.5 * g.areaNuE0TE0T{nn, np, m}.', RtildeT).';
        else
          ret{m} = ret{m} + kronVec(0.5 * bsxfun(@times, g.markE0TE0T{nn, np}, areaNuMarkE0T{m}(:, nn)).', RtildeT).';
        end % if
      end % for np
      % Diagonal blocks
      for l = 1 : dataN
        ret{m} = ret{m} + kron(spdiags(0.5 * areaNuMarkE0T{m}(:, nn) .* dataDisc{m}(:, l), 0, K, K), refEdgePhiIntPhiIntPhiInt(:, :, l, nn));
      end % for l
    end % for m
  end % for nn
  
end % if
end % function


function ret = assembleMatEdgePhiPhiFuncDiscScalarNu(g, refEdgePhiIntPhiIntPhiInt, refEdgePhiIntPhiExtPhiExt, dataDisc, areaNuMarkE0T)
[K, dataN] = size(dataDisc);  N = size(refEdgePhiIntPhiIntPhiInt, 1);  nEdges = size(g.E0T, 2);

ret = { sparse(K*N, K*N), sparse(K*N, K*N) };

if numel(g.markE0TE0T) == nEdges % mapping from nn to np implicitly given
  
  validateattributes(refEdgePhiIntPhiExtPhiExt, {'numeric'}, {'size', [N N N nEdges]}, mfilename, 'refEdgePhiIntPhiExtPhiExt');
  
%   for n = 1 : nEdges
%     % Coefficient-function for off-diagonal blocks
%     RtildeT = zeros(K*N, N);
%     for l = 1 : dataN
%       RtildeT = RtildeT + kron(dataDisc(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, n).');
%     end % for l
%     for m = 1 : 2
%       % Off-diagonal blocks
%       ret{m} = ret{m} + kronVec((spdiags(0.5 * areaNuMarkE0T{m}(:, n), 0, K, K) * g.markE0TE0T{n}).', RtildeT).';
%       % Diagonal blocks
%       for l = 1 : dataN
%         ret{m} = ret{m} + kron(spdiags(0.5 * areaNuMarkE0T{m}(:, n) .* dataDisc(:, l), 0, K, K), ...
%                                refEdgePhiIntPhiIntPhiInt(:, :, l, n));
%       end % for l
%     end % for m
%   end  % for n
  
  % Diagonal blocks
  for m = 1 : 2
    RtildeT = zeros(K*N, N);
    for n = 1 : nEdges
      for l = 1 : dataN
        RtildeT = RtildeT + kron(0.5 * areaNuMarkE0T{m}(:, n) .* dataDisc(:, l), refEdgePhiIntPhiIntPhiInt(:, :, l, n));
      end % for l
    end % for m
    ret{m} = kronVec(speye(K, K), RtildeT);
  end  % for n
  
  % Off-diagonal blocks
  for n = 1 : nEdges
    % Coefficient-function for off-diagonal blocks
    RtildeT = zeros(K*N, N);
    for l = 1 : dataN
      RtildeT = RtildeT + kron(dataDisc(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, n).');
    end % for l
    for m = 1 : 2
      ret{m} = ret{m} + kronVec((spdiags(0.5 * areaNuMarkE0T{m}(:, n), 0, K, K) * g.markE0TE0T{n}).', RtildeT).';
    end % for m
  end  % for n
    
  
else % mapping from nn to np explicitly given
  
  validateattributes(refEdgePhiIntPhiExtPhiExt, {'numeric'}, {'size', [N N N nEdges nEdges]}, mfilename, 'refEdgePhiIntPhiExtPhiExt');

  for nn = 1 : nEdges
    % Off-diagonal blocks
    for np = 1 : nEdges
      RtildeT = zeros(K*N, N);
      for l = 1 : dataN
        RtildeT = RtildeT + kron(dataDisc(:, l), refEdgePhiIntPhiExtPhiExt(:, :, l, nn, np).');
      end % for
      if isfield(g, 'areaNuE0TE0T')
        ret{1} = ret{1} + kronVec(0.5 * g.areaNuE0TE0T{nn, np, 1}.', RtildeT).';
        ret{2} = ret{2} + kronVec(0.5 * g.areaNuE0TE0T{nn, np, 2}.', RtildeT).';
      else
        ret{1} = ret{1} + kronVec(0.5 * bsxfun(@times, g.markE0TE0T{nn, np}, areaNuMarkE0T{1}(:, nn)).', RtildeT).';
        ret{2} = ret{2} + kronVec(0.5 * bsxfun(@times, g.markE0TE0T{nn, np}, areaNuMarkE0T{2}(:, nn)).', RtildeT).';
      end % if
    end % for
    % Diagonal blocks
    for l = 1 : dataN
      ret{1} = ret{1} + kron(spdiags(0.5 * areaNuMarkE0T{1}(:, nn) .* dataDisc(:,l), 0, K, K), refEdgePhiIntPhiIntPhiInt(:, :, l, nn));
      ret{2} = ret{2} + kron(spdiags(0.5 * areaNuMarkE0T{2}(:, nn) .* dataDisc(:,l), 0, K, K), refEdgePhiIntPhiIntPhiInt(:, :, l, nn));
    end % for
  end % for
  
end % if
end % function
