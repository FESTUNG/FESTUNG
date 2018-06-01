% It computes a rotation matrix, according to a specified variant.
%
%===============================================================================
%> @file
%>
%> @brief It computes a rotation matrix, according to a specified variant.
%===============================================================================
%>
%> @brief It computes a rotation matrix, according to a specified variant.
%> 
%> This function is called by the vector vertex based limiters to compute 
%> the limiting directions (cf. <code>applyVectorVertexBasedLimiter()<code>
%> and <code>computeSequentialVectorVertexBasedLimiter()</code>).
%>
%> @param  g            The lists describing the geometric and topological 
%>                      properties of a triangulation (see 
%>                      <code>generateGridData()</code>) 
%>                      @f$[1 \times 1 \text{ struct}]@f$
%> @param  gradient     The gradient of the unlimited vector field 
%>                      @f$\mathbf{u}_h@f$. @f$[2 \times 2 \mathrm{cell}]@f$
%> @param  variant      The variant that is used to compute rotation
%>                      matrices for directional vector limiting.
%> @param  angle        The rotation angle to be used in case a fixed
%>                      transformation is desired to compute the limiting directions.
%> @retval ret          The computed rotation matrixces for each element.
%>                      @f$[2 \times 2 \mathrm{cell}]@f$
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
function ret = computeRotationMatrix(g, gradient, variant, angle)
switch variant
  case 'fixed'
    K = g.numT;
    c = cos(angle); s = sin(angle);
    ret = {c*ones(K,1), s*ones(K,1); -s*ones(K,1), c*ones(K,1)}; % backtransformation to unrotated grid
  case 'transformation'
    ret = computeLeftSVDRotationMatrix( {g.B(:,1,1), g.B(:,1,2); g.B(:,2,1), g.B(:,2,2)} );
  case 'gradient'
    ret = computeLeftSVDRotationMatrix(gradient);
  case 'Gram-Schmidt'
    ret = computeBasisFromGramSchmidt({g.B(:,1,1), g.B(:,2,1)}, {g.B(:,1,2), g.B(:,2,2)});
  case 'kernel'
    ret = computeKernelBasedRotationMatrix(gradient);
  otherwise
    error('Unknown variant.')
end % switch
end % function
