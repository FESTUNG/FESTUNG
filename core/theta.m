% Mapping from edge nn to edge np in the reference triangle.
%
%===============================================================================
%> @file
%>
%> @brief Mapping from edge @f$\hat{E}_{n^-}@f$ to edge @f$\hat{E}_{n^+}@f$ in
%>        the reference triangles.
%===============================================================================
%>
%> @brief Computes the reference coordinates
%>        @f$\hat{\mathbf{x}}_p = (\hat{x}^1,\hat{x}^2)^T@f$ as given by the 
%>        mapping @f$\hat{\mathbf{x}}_p = 
%>        \hat{\mathbf{\vartheta}}_{n^-n^+}(\hat{\mathbf{x}})@f$.
%>
%> The mapping 
%> @f$\hat{\mathbf{\vartheta}}_{n^-n^+} : \hat{E}_{n^-} \mapsto \hat{E}_{n^+}@f$ 
%> maps coordinates from edge @f$\hat{E}_{n^-}, n^- \in \{1,2,3\}@f$ to edge 
%> @f$\hat{E}_{n^+}, n^+ \in \{1,2,3\}@f$ in the reference triangle 
%> @f$\hat{T} = \{ (0,0), (1,0), (0,1) \}@f$.
%>
%> For an arbitrary index pair @f$k^-,k^+@f$ it is defined as
%> @f[
%> \hat{\mathbf{\vartheta}}_{n^-n^+}(\hat{\mathbf{x}}) = 
%>   \mathbf{F}^{-1}_{k^+} \circ \mathbf{F}_{k^-}(\hat{\mathbf{x}}) \,,
%> @f]
%> with @f$\mathbf{F}_k : \hat{T} \mapsto T_k@f$ the affine mapping from
%> reference triangle @f$\hat{T}@f$ to physical triangle @f$T@f$, defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,.
%> @f]
%> This can be boiled down to nine possible maps between the sides of the
%> reference triangle, the closed form expressions of which read as
%>
%> @f[
%>  \hat{\mathbf{\vartheta}}_{11}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} 1 - \hat{x}^1 \\ 1 - \hat{x}^2 \end{bmatrix} \,,
%>  \hat{\mathbf{\vartheta}}_{21}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} 1 - \hat{x}^2 \\ \hat{x}^2 \end{bmatrix} \,, 
%>  \hat{\mathbf{\vartheta}}_{31}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} \hat{x}^1 \\ 1 - \hat{x}^1 \end{bmatrix} \,, 
%> @f]
%> @f[
%>  \hat{\mathbf{\vartheta}}_{12}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} 0 \\ \hat{x}^2 \end{bmatrix} \,,
%>  \hat{\mathbf{\vartheta}}_{22}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} 0 \\ 1 - \hat{x}^2 \end{bmatrix} \,, 
%>  \hat{\mathbf{\vartheta}}_{32}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} 0 \\ \hat{x}^1 \end{bmatrix} \,, 
%> @f]
%> @f[
%>  \hat{\mathbf{\vartheta}}_{13}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} \hat{x}^1 \\ 0 \end{bmatrix} \,,
%>  \hat{\mathbf{\vartheta}}_{23}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} \hat{x}^2 \\ 0 \end{bmatrix} \,, 
%>  \hat{\mathbf{\vartheta}}_{33}:
%>    \begin{bmatrix} \hat{x}^1 \\ \hat{x}^2 \end{bmatrix} \mapsto
%>    \begin{bmatrix} 1 - \hat{x}^1 \\ 0 \end{bmatrix} \,.
%> @f]
%> All maps @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ reverse the edge orientation
%> because an edge shared by triangles @f$T^-@f$ and @f$T^+@f$ will always have
%> different orientations when mapped by @f$\mathbf{F}_{k^-}@f$ and 
%> @f$\mathbf{F}_{k^+}@f$; this occurs due to the counter-clockwise vertex
%> orientation consistently maintained throughout the mesh.
%>
%> @param  nn   The index @f$n^-@f$ of edge @f$\hat{E}_{n^-}@f$ in the 
%>              reference triangle @f$\hat{T}@f$ @f$[\text{scalar}]@f$
%> @param  np   The index @f$n^+@f$ of edge @f$\hat{E}_{n^+}@f$ in the 
%>              reference triangle @f$\hat{T}@f$ @f$[\text{scalar}]@f$
%> @param  X1,X2 The coordinates @f$\hat{\mathbf{x}}=(\hat{x}^1,\hat{x}^2)^T@f$
%>               on edge @f$\hat{E}_{n^-}@f$.
%> @retval XP1,XP2 The coordinates 
%>               @f$\hat{\mathbf{x}}_p=(\hat{x}_p^1,\hat{x}_p^2)^T@f$
%>               on edge @f$\hat{E}_{n^+}@f$.
%>               It holds <code>size(XP1) == size(XP2) == size(X1) == size(X2)</code>.
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
function [XP1, XP2] = theta(nn, np, X1, X2)
assert(1 <= nn && 3 >= nn, 'nn must be in [1,3]')
assert(1 <= np && 3 >= np, 'np must be in [1,3]')
assert(isequal(size(X1), size(X2)), 'X1 and X2 must have the same size')
switch nn
  case 1
    switch np
      case 1,  XP1 = 1-X1;  XP2 = 1-X2;
      case 2,  XP1 = zeros(size(X1));  XP2 = X2;
      case 3,  XP1 = X1;  XP2 = zeros(size(X1));
    end % switch
  case 2
    switch np
      case 1,  XP1 = 1-X2;  XP2 = X2;
      case 2,  XP1 = zeros(size(X1));  XP2 = 1-X2;
      case 3,  XP1 = X2;  XP2 = zeros(size(X1));
    end % switch
  case 3
    switch np
      case 1,  XP1 = X1;  XP2 = 1-X1;
      case 2,  XP1 = zeros(size(X1));  XP2 = X1;
      case 3,  XP1 = 1-X1;  XP2 = zeros(size(X1));
    end % switch
end % switch
end % function
