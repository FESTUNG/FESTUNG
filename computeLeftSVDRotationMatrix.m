% Computes the first rotation matrix of the singular value decomposition of a 
% 2x2 matrix using an algebraic formula based on the polar decomposition.

%===============================================================================
%> @file
%>
%> @brief Computes the first rotation matrix of the singular value decomposition
%>        of a 2x2 matrix using an algebraic formula based on the polar 
%>        decomposition.
%===============================================================================
%>
%> @brief Computes the first rotation matrix of the singular value decomposition
%>        of a 2x2 matrix using an algebraic formula based on the polar 
%>        decomposition.
%>
%> It uses an analytical formula to compute the left rotation matrices in 
%> the singular value decomposition of a list of given 2x2 matrices.
%> Since such matrices are not unique (not even up to the factor -1) in the
%> case of equal or zero singular values, the robustness of this method can
%> be optimized by calling it from 
%> <code>computeKernelBasedRotationMatrix()</code>. 
%> In the case of equal singular values, the rotation matrix is 
%> automatically chosen such that the product of input and output matrices
%> is a mutiple of unity. 
%> Note that the influence of the chosen tolerance value may be critical.
%> The auxiliary routine
%> <code>matrixMultiplication()</code>
%> is utilized to compute products of lists of 2x2 matrices.
%>
%> @param  mat        A cell describing a list of 2x2 of matrices.
%>                    @f$[2 \times 2 \text{ cell}]@f$
%> @retval ret        A cell describing the list of computed rotation 
%>                    matrices.
%>                    @f$[2 \times 2 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2018.
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
function ret = computeLeftSVDRotationMatrix(mat)

validateattributes(mat, {'cell'}, {'size' [2 2]}, mfilename, 'mat')

K = length(mat{1,1});
tol = 1e-12;

ret = cell(2,2);
ret{1,1} =  ones(K,1); ret{1,2} = zeros(K,1);
ret{2,1} = zeros(K,1); ret{2,2} =  ones(K,1);


E = 0.5 * (mat{1,1} + mat{2,2});
F = 0.5 * (mat{1,1} - mat{2,2});
G = 0.5 * (mat{2,1} + mat{1,2});
H = 0.5 * (mat{2,1} - mat{1,2});
Q = sqrt(E.^2 + H.^2);
R = sqrt(F.^2 + G.^2);
sx = Q+R;
sy = Q-R;

ind1 = sx > tol;
ind2 = abs(sy) > tol;
ind3 = abs(R) > tol;

% ind = ~ind1 & ~ind2;
% use unity, i.e. nothing to do

ind = ind1 & ~ind2; % sx > sy = 0
if any(ind)
  vec = [mat{1,1}(ind), mat{2,1}(ind)];
  norm = sqrt(vec(:,1).^2 + vec(:,2).^2);
  vec = [vec(:,1) ./ norm, vec(:,2) ./ norm];
  ret{1,1}(ind) =  vec(:,1);
  ret{1,2}(ind) =  vec(:,2);
  ret{2,1}(ind) = -vec(:,2);
  ret{2,2}(ind) =  vec(:,1);
end % if

ind = ind1 & ~ind3; % sx = sy
if any(ind)
  coef = 1 ./ sx(ind);
  ret{1,1}(ind) = coef .* mat{1,1}(ind);
  ret{1,2}(ind) = coef .* mat{2,1}(ind);
  ret{2,1}(ind) = coef .* mat{1,2}(ind);
  ret{2,2}(ind) = coef .* mat{2,2}(ind);
end % if

ind = ind1 & ind2 & ind3; % 0 < sx ~= sy > 0
if any(ind)
  a1 = atan2(G(ind), F(ind));
  a2 = atan2(H(ind), E(ind));
  theta = 0.5 * (a2 - a1);
  phi   = 0.5 * (a2 + a1);
  S = 2*(sy(ind) >= 0)  - 1;

  c = cos(phi);
  s = sin(phi);
  U = {c, s; -S.*s, S.*c};

  c = cos(theta);
  s = sin(theta);
  V = {c, s; -s, c};

  D = matrixMultiplication(matrixMultiplication(U, {mat{1,1}(ind), mat{1,2}(ind); mat{2,1}(ind), mat{2,2}(ind)}), V);

  if max(max(abs(D{1,2}), abs(D{2,1}))) > tol
    error('SVD results in non-diagonal matrix D.');
  end % if
  if max(max( abs(D{1,1} - sx(ind)), abs(D{2,2} - abs(sy(ind))) )) > tol
    error('SVD results in wrong diagonal matrix D.');
  elseif min(sx(ind)) < -tol
    error('SVD results in non positive diagonal matrix D.');
  end % if

  ret{1,1}(ind) = U{1,1};
  ret{1,2}(ind) = U{1,2};
  ret{2,1}(ind) = U{2,1};
  ret{2,2}(ind) = U{2,2};
end % if

end % function

function ret = matrixMultiplication(A, B)

ret = { A{1,1}.*B{1,1}+A{1,2}.*B{2,1}, A{1,1}.*B{1,2}+A{1,2}.*B{2,2}; ...
				A{2,1}.*B{1,1}+A{2,2}.*B{2,1}, A{2,1}.*B{1,2}+A{2,2}.*B{2,2} };

end % function
