% Computes a modified Kronecker product with a vector of matrices as second
% operand.

%===============================================================================
%> @file
%>
%> @brief Computes a modified Kronecker product with a vector of matrices 
%>        as second operand.
%===============================================================================
%>
%> @brief Computes a modified Kronecker product with a vector of matrices 
%>        as second operand.
%> 
%> It evaluates the operation 
%> @f$\mathsf{K} = \mathsf{A} \otimes_\mathrm{V} \mathsf{B}@f$, with
%> @f$\mathsf{A} = [a_{ij}] \in \mathbb{R}^{m_a \times n_a}@f$, 
%> @f$\mathsf{B} = [b_{kl}] \in \mathbb{R}^{m_b \times n_b}@f$,
%> @f$m_b = r m_a, r\in\mathbb{N}@f$ and the operation defined as
%> @f[
%>  \mathsf{A} \otimes_\mathrm{V} \mathsf{B} \;:=\;
%>  \left[a_{ij} \left[\mathsf{B}\right]_{(i-1)r\,:\,ir, :} \right]
%>  = \mathsf{K} \in \mathbb{R}^{m_b \times n_a n_b}\;.
%> @f]
%> It can be understood as a Kronecker product that takes a different right
%> hand side for every row of the left hand side.  The input matrices
%> @f$\mathsf{A}@f$ and @f$\mathsf{B}@f$ may be sparse or dense.  If at
%> least one input argument is sparse, @f$\mathsf{K}@f$ will be sparse.
%>
%> @par Example
%> @parblock
%> @f[
%>  \mathsf{A} := \begin{bmatrix}1 & 2 \\ 3 & 4\end{bmatrix}~,\quad
%>  \mathsf{B} := \begin{bmatrix}1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 1 & 1 \\ 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}~.
%> @f]
%>
%> Here, @f$r = 3@f$.  It follows
%> 
%> @f[
%>  \mathsf{K} = \begin{bmatrix}1 & 1 & 1 & 2 & 2 & 2 \\ 
%>                               1 & 1 & 1 & 2 & 2 & 2 \\
%>                               1 & 1 & 1 & 2 & 2 & 2 \\
%>                               3 & 0 & 0 & 0 & 0 & 0 \\
%>                               0 & 3 & 0 & 0 & 0 & 0 \\
%>                               0 & 0 & 3 & 0 & 0 & 0
%>                \end{bmatrix}~. 
%> @f]
%> @endparblock
%>
%> @note The number of rows of the second operand @f$\mathsf{B}@f$ must be
%>       the same or a multiple of the number of rows of the first operand
%>       @f$A@f$.
%>
<<<<<<< HEAD
=======
%> 
>>>>>>> github_new
%> @param A  The first operand. @f$[m_a \times n_a]@f$
%> @param B  The second operand. @f$[m_b \times n_b]@f$
%> @retval K The result of the operation. @f$[m_b \times n_a n_b]@f$
%>
%> This file is part of FESTUNG
%>
<<<<<<< HEAD
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
=======
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
>>>>>>> github_new
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
function K = kronVec(A, B)
if ~ismatrix(A) || ~ismatrix(B)
  error(message('MATLAB:kron:TwoDInput'));
end
[ma,na] = size(A);
[mb,nb] = size(B);
if ~mod(mb, ma) == 0
  error('The number of rows in B must be a multiple of the number of rows in A!');
end
mc = mb / ma;
if ~issparse(A) && ~issparse(B)
  % Both inputs full, result is full.
  A = reshape(A, [1 ma 1 na]);
  B = reshape(B, [mc ma nb 1]);
  K = reshape(bsxfun(@times, A, B), [mb na*nb]);
else
  % Replicate rows of A
  [i2, j2, v2] = find(kron(A, ones(mc, 1)));
  % Compute row and column vectors
  ik = repmat(i2, [1 nb]);
  jk = bsxfun(@plus, nb * (j2 - 1), 1 : nb);
  % Select and replicate blocks from B and multiply with A
  sk = bsxfun(@times, v2, B(i2, :));
  % Build result matrix
  K = sparse(ik, jk, sk, mb, na * nb);
end % if
end % function
