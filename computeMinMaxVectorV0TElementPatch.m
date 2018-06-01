% Assembles two matrices with minimum or maximum values, respectively. 
% The computed values are scalar products of a limiting direction and 
% differences of centroid values of the patch of elements surrounding each 
% vertex of each element of a vector field. 

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices with minimum or maximum values, respectively. 
%>        The computed values are scalar products of a limiting direction and 
%>        differences of centroid values of the patch of elements surrounding
%>        each vertex of each element of a vector field. 
%===============================================================================
%>
%> @brief Assembles two matrices with minimum or maximum values, respectively. 
%>        The computed values are scalar products of a limiting direction and 
%>        differences of centroid values of the patch of elements surrounding
%>        each vertex of each element of a vector field. 
%>
%> For each vertex @f$\mathbf{x}_{ki}, i\in\{1,2,3\}@f$ of each triangle
%> @f$T_k@f$ it computes the minimum and maximum centroid values of
%> @f$\mathbf{q}_j^k \cdot (\mathbf{u}_{l\mathrm{c}} - \mathbf{u}_{k\mathrm{c}})@f$
%> for a vector field @f$\mathbf{u}_h(\mathbf{x})@f$ and a limiting
%> direction @f$\mathbf{q}_j^k@f$ for each element.
%> @f[
%> \forall T_k \in \mathcal{T}_h \,, \forall i\in\{1,2,3\} \,, \quad
%> \mathbf{u}_{ki}^\mathrm{min} := \min_{\{ T_l \in \mathcal{T}_h \,|\, \mathbf{x}_{ki} \in T_l \}}
%> mathbf{q}_j^k \cdot (\mathbf{u}_{l\mathrm{c}} - \mathbf{u}_{k\mathrm{c}}) \,,\quad
%> \mathbf{u}_{ki}^\mathrm{max} := \max_{\{ T_l \in \mathcal{T}_h \,|\, \mathbf{x}_{ki} \in T_l \}}
%> mathbf{q}_j^k \cdot (\mathbf{u}_{l\mathrm{c}} - \mathbf{u}_{k\mathrm{c}}) \,,\quad
%> @f]
%> where @f$\mathbf{u}_{l\mathrm{c}} = \mathbf{u}_h(\mathbf{x}_{lc})@f$ is 
%> the function value of the vector field in the centroid @f$\mathbf{x}_{lc}@f$ of triangle 
%> @f$T_l@f$.
%>
%> @param  g           The lists describing the geometric and topological 
%>                     properties of a triangulation (see 
%>                     <code>generateGridData()</code>)
%>                     @f$[1 \times 1 \text{ struct}]@f$
%> @param  vector      A list of directions used in scalar products with
%>                     the vector field that is being limited.
%> @param  valCentroid The centroid values of the vector field to be limited.
%>                     @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  markV0TbdrD <code>logical</code> array that marks each triangles
%>                     (Dirichlet boundary) vertices on which the values given
%>                     in dataV0T should be included in the computation of
%>                     minimum/maximum @f$[K \times 3]@f$
%> @param  dataV0T     A list of values to be included in the determination
%>                     of minimum/maximum. @f$[2 \times 1 \mathrm{cell}]@f$
%> @retval ret         The assembled matrices @f$[K \times 3]@f$ as
%>                     @f$[2 \times 1 \mathrm{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2018
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
function minMaxV0T = computeMinMaxVectorV0TElementPatch(g, vector, valCentroid, markV0TbdrD, dataV0T)

% Check function arguments that are directly used
validateattributes(vector, {'cell'}, {'size', [1 2]}, mfilename, 'vector');
validateattributes(valCentroid, {'cell'}, {'size', [2 1]}, mfilename, 'valCentroid');
validateattributes(markV0TbdrD, {'logical'}, {'size', [g.numT 3]}, mfilename, 'markV0TbdrD');
validateattributes(dataV0T, {'cell'}, {'size', [2 1]}, mfilename, 'dataV0T');

% Include boundary value on Dirichlet vertices into limiting procedure
valD{1} = NaN(g.numT, 3);  valD{2} = NaN(g.numT, 3);
valD{1}(markV0TbdrD) = dataV0T{1}(markV0TbdrD);  valD{2}(markV0TbdrD) = dataV0T{2}(markV0TbdrD);
valD = repmat(vector{1}, [1 3]) .* valD{1} + repmat(vector{2}, [1 3]) .* valD{2};

if isfield(g, 'markV0TT0V')
  minMaxV0T = computeMinMaxVectorV0T_withMarkV0TT0V(g, vector, valCentroid, valD);
else
  minMaxV0T = computeMinMaxVectorV0T_noMarkV0TT0V(g, vector, valCentroid, valD);
end
end % function
%
%===============================================================================
%> @brief Helper function for the case that computeMinMaxVectorV0TElementPatch()
%> was called with a precomputed field markV0TT0V in parameter g.
%
function minMaxV0T = computeMinMaxVectorV0T_withMarkV0TT0V(g, vector, valCentroid, valD)
% Initialize return value: a 2x1-cell with a Kx3-array of minimum and 
% maximum values per vertex of each triangle
minMaxV0T = cell(2,1); 
minMaxV0T{1} = zeros(g.numT, 3); minMaxV0T{2} = zeros(g.numT, 3);

for i = 1 : 3
  diff = bsxfun(@times, bsxfun(@times, g.markV0TT0V{i}, vector{1}), valCentroid{1}') + bsxfun(@times, bsxfun(@times, g.markV0TT0V{i}, vector{2}), valCentroid{2}') - bsxfun(@times, g.markV0TT0V{i}, valCentroid{1} .* vector{1} + valCentroid{2} .* vector{2});
  shiftValToPos = abs(min(diff, [], 2)) + 1;
  shiftValToNeg = abs(max(diff, [], 2)) + 1;

  % Compute min/max per vertex as min/max of centroid values per vertex
  minMaxV0T{1}(:, i) = min( min( diff - bsxfun(@times, g.markV0TT0V{i}, shiftValToNeg), [], 2 ) + shiftValToNeg, valD(:,i) );
  minMaxV0T{2}(:, i) = max( max( diff + bsxfun(@times, g.markV0TT0V{i}, shiftValToPos), [], 2 ) - shiftValToPos, valD(:,i) );
end %for
end % function
%
%===============================================================================
%> @brief Helper function for the case that computeMinMaxVectorV0TElementPatch()
%> was called with no precomputed field markV0TT0V in parameter g.
%
function minMaxV0T = computeMinMaxVectorV0T_noMarkV0TT0V(g, vector, valCentroid, valD)
% Initialize return value: a 2x1-cell with a Kx3-array of minimum and 
% maximum values per vertex of each triangle
minMaxV0T = cell(2,1); 
minMaxV0T{1} = zeros(g.numT, 3); minMaxV0T{2} = zeros(g.numT, 3);

for i = 1 : 3
  % Mark all elements sharing i-th vertex
  markV0TT0V = g.markV0TV0T{i, 1} | g.markV0TV0T{i, 2} | g.markV0TV0T{i, 3}; 

  % Fix a bug in GNU Octave 4.0.0's implementation of sparse matrix concatenation
  if exist('OCTAVE_VERSION','builtin')
    markV0TT0V = markV0TT0V + 0 * speye(size(markV0TT0V, 1), size(markV0TT0V, 2));
  end % if
  
  diff = bsxfun(@times, bsxfun(@times, markV0TT0V, vector{1}), valCentroid{1}') + bsxfun(@times, bsxfun(@times, markV0TT0V, vector{2}), valCentroid{2}') - bsxfun(@times, markV0TT0V, valCentroid{1} .* vector{1} + valCentroid{2} .* vector{2});
  shiftValToPos = abs(min(diff, [], 2)) + 1;
  shiftValToNeg = abs(max(diff, [], 2)) + 1;

  % Compute min/max per vertex as min/max of centroid values per vertex
  minMaxV0T{1}(:, i) = min( min( diff - bsxfun(@times, markV0TT0V, shiftValToNeg), [], 2 ) + shiftValToNeg, valD(:,i) );
  minMaxV0T{2}(:, i) = max( max( diff + bsxfun(@times, markV0TT0V, shiftValToPos), [], 2 ) - shiftValToPos, valD(:,i) );
end % for
end % function
