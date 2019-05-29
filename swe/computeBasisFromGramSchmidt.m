% Performs the Gram-Schmidt process for a list of pairs of vectors with two
% entries and puts the resulting vectors row-wise into a 2x2 matrix.

%===============================================================================
%> @file
%>
%> @brief Performs the Gram-Schmidt orthogonalization process for a list of
%>        pairs of vectors with two entries and puts the resulting vectors
%>        row-wise into a 2x2 matrix.
%===============================================================================
%>
%> @brief Performs the Gram-Schmidt orthogonalization process for a list of
%>        pairs of vectors with two entries and puts the resulting vectors
%>        row-wise into a 2x2 matrix.
%>
%> First, the following distinction is made:
%> If both vectors in a corresponding pair are zero, the output is unity.
%> If either of the pair of input vectors is zero but the other is not, 
%> the latter is normalized and its orthogonal complement is used as second
%> basis vector.
%> Otherwise the regular Gram-Schmidt process is performed.
%> Note that the influence of the chosen tolerance value may be critical.
%>
%> @param  v1        A cell describing a list of of vectors in \mathbb{R^2}.
%>                   @f$[2 \times 1 \text{ cell}]@f$
%> @param  v2        A cell describing a list of of vectors in \mathbb{R^2}.
%>                   @f$[2 \times 1 \text{ cell}]@f$
%> @retval ret       A 2x2 cell describing the list of computed rotation 
%>                   matrices, formed row-wisefrom the orthonormal basis
%>                   computed in the Gram-Schmidt process.
%>                   @f$[2 \times 2 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = computeBasisFromGramSchmidt(v1, v2)

tol = 1e-8;

K = length(v1{1});

ret = cell(2,2);
ret{1,1} =  ones(K,1); ret{1,2} = zeros(K,1);
ret{2,1} = zeros(K,1); ret{2,2} =  ones(K,1);

ind1 = (abs(v1{1}) < tol) & (abs(v1{2}) < tol);
ind2 = (abs(v2{1}) < tol) & (abs(v2{2}) < tol);

% ind = ind1 & ind2; nothing to do for zero matrix

ind = ind1 & ~ind2;
if any(ind)
  
  norm = sqrt(v2{1}(ind).^2 + v2{2}(ind).^2);
  w1 = { v2{1}(ind) ./ norm; v2{2}(ind) ./ norm };

  ret{1,1}(ind) = -w1{2}; ret{1,2}(ind) = w1{1};
  ret{2,1}(ind) =  w1{1}; ret{2,2}(ind) = w1{2};
  
end % if

ind = ~ind1 & ind2;

if any(ind)
  norm = sqrt(v1{1}(ind).^2 + v1{2}(ind).^2);
  w1 = { v1{1}(ind) ./ norm; v1{2}(ind) ./ norm };
  
  ret{1,1}(ind) =  w1{1}; ret{1,2}(ind) = w1{2};
  ret{2,1}(ind) = -w1{2}; ret{2,2}(ind) = w1{1};
  
end % if

ind = ~ind1 & ~ind2;
if any(ind)
  norm = sqrt(v1{1}(ind).^2 + v1{2}(ind).^2);
  w1 = { v1{1}(ind) ./ norm; v1{2}(ind) ./ norm };
  
  prod = v2{1}(ind) .* w1{1} + v2{2}(ind) .* w1{2};
  w2 = { v2{1}(ind) - prod .* w1{1}; v2{2}(ind) - prod .* w1{2} };
  
  norm = sqrt(w2{1}.^2 + w2{2}.^2);
  w2 = { w2{1} ./ norm; w2{2} ./ norm };
  
  ret{1,1}(ind) = w1{1}; ret{1,2}(ind) = w1{2};
  ret{2,1}(ind) = w2{1}; ret{2,2}(ind) = w2{2};
  
end % if

end % function
