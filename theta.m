% Mapping from edge nn to edge np in the reference triangle.
%
%> @file theta.m
%>
%> @brief Mapping from edge @f$\hat{E}_{n^-}@f$ to edge @f$\hat{E}_{n^+}@f$ in
%>        the reference triangles
%>
%> <code>[XP1,XP2] = theta(nn, np, X1, X2)</code> returns the reference coordinates
%> @f$\hat{\mathbf{x}}_p = (\hat{x}^1,\hat{x}^2)^T@f$ as given by the mapping
%> @f$\hat{\mathbf{x}}_p=\hat{\mathbf{\vartheta}}_{n^-n^+}(\hat{\mathbf{x}})@f$.
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
%>
%> @param  nn   The index @f$n^-@f$ of edge @f$\hat{E}_{n^-}@f$ in the 
%>              reference triangle @f$\hat{T}@f$ @f$[\text{scalar}]@f$
%> @param  np   The index @f$n^+@f$ of edge @f$\hat{E}_{n^+}@f$ in the 
%>              reference triangle @f$\hat{T}@f$ @f$[\text{scalar}]@f$
%> @param  S  The parameter @f$s@f$ of the mapping. Can be a vector, e.g.,
%>            to compute the mapping for multiple values in one function call.
%> @retval X1,X2  reference coordinates
%>                @f$\hat{\mathbf{x}} = (\hat{x}^1,\hat{x}^2)^T@f$.
%>                It holds <code>size(X1) == size(X2) == size(S)</code>.
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
