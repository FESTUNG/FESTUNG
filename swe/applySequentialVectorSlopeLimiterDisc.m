% Applies a sequential slope limiter to a discrete function and a two 
% component vector of products of functions with the former, all given in 
% (modal) DG basis.

%===============================================================================
%> @file
%>
%> @brief Applies a sequential slope limiter to a discrete function and a 
%>        two component vector of products of functions with the former, 
%>        all given in (modal) DG basis.
%===============================================================================
%>
%> @brief Applies a sequential slope limiter to a discrete function and a 
%>        two component vector of products of functions with the former, 
%>        all given in (modal) DG basis.
%>
%> It applies the sequential slope limiting operator for a two component 
%> vector field
%> @f[
%>  \mathsf{\Phi} := \mathsf{M}^{-1} \mathsf{M}^\mathrm{DG,Taylor} 
%>  \mathsf{\Phi}^\mathrm{Taylor,seq} \left(\mathsf{M}^\mathrm{DG,Taylor}\right)^{-1} \mathsf{M}
%> @f]
%> to a discrete function @f$\xi_h: \Omega \rightarrow \mathbb{R}@f$ given as
%> representation matrix, e.g., obtained by L2-projection of a continuous
%> function (cf. <code>projectFunctCont2DataDisc()</code>),
%> and a vector field @f$\mathbf{u}_h: \Omega \rightarrow \mathbb{R}^d@f$
%> obtained from the quotient @f$\frac{\mathbf{u}_h}{\xi_h-b_h}@f$, where 
%> @f$b_h@f$ is a discrete function given by its DG coefficients.
%> @f$\mathsf{M}@f$ is the mass matrix in DG basis (cf.
%> <code>assembleMatElemPhiPhi()</code>),
%> @f$\mathsf{M}^\mathrm{DG,Taylor}@f$ is the transformation matrix from
%> (modal) DG to Taylor basis (cf.
%> <code>assembleMatElemPhiDGPhiTaylor()</code>), and
%> @f$\mathsf{\Phi}^\mathrm{Taylor,seq}@f$ is the sequential slope limiting
%> operator in Taylor basis 
%> (cf. <code>computeSequentialVertexBasedLimiter()</code>) for a two 
%> component vector field.
%>
%> @param  g                The lists describing the geometric and topological 
%>                          properties of a triangulation (see 
%>                          <code>generateGridData()</code>) 
%>                          @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc         The representation matrix of the unlimited functions 
%>                          @f$\xi_h@f$ and @f$\mathbf{uH}_h@f$.
%> @param  zbDisc           The representation matrix of the problem
%>                          bathymetry @f$b_h@f$ @f$[K \times N]@f$.
%> @param  markV0TbdrD      <code>cell</code> that for each unknown marks
%>                          each triangles (Dirichlet boundary) vertices on 
%>                          which additional function values should be
%>                          considered during the  slope limiting routine. 
%>                          @f$[3 \times 1]@f$
%> @param  dataV0T          The function values for each unknown for 
%>                          (Dirichlet boundary) vertices specified by 
%>                          <code>markV0TbdrD</code> @f$[3 \times 1]@f$ <code>cell</code>.
%> @param  globM            The mass matrix @f$\mathsf{M}@f$ in DG basis.
%>                          @f$[KN \times KN]@f$.
%> @param  globMDiscTaylor  The transformation matrix @f$\mathsf{M}^\mathrm{DG,Taylor}@f$.
%>                          @f$[KN \times KN]@f$.
%> @param  basesOnQuad      A struct containing precomputed values of (Taylor) basis
%>                          functions on quadrature points. Must provide at
%>                          least phiTaylorV0T.
%> @param  type             The type of slope limiter to be used. [<code>string</code>]
%> @param  variant          The variant that is used to compute rotation
%>                          matrices for directional vector limiting.
%> @param  angle            The rotation angle to be used in case a fixed
%>                          transformation is desired to compute the limiting directions.
%> @retval dataDisc         The representation matrix of the limited functions
%>                          @f$\mathsf{\Phi}c_h@f$. @f$[K \times N \times d+1]@f$
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
function dataDisc = applySequentialVectorSlopeLimiterDisc(g, dataDisc, zbDisc, markV0TbdrD, dataV0T, globM, globMDiscTaylor, basesOnQuad, type, variant, angle)

% expects d+1 (d > 0) components in cells dataDisc, markV0TbdrD, dataV0T
% each component of markV0TbdrD, dataV0T may be empty or meaningful for
% as long as it is consistent with each other
% Limits dataDisc{1} and then dataDisc{1+i} according to the maximum
% principle for the quotinet of the variabes corresponding to components
% 1+i and 1.

numProdVar = length(dataDisc);
dataTaylor = cell(numProdVar, 1);
% Project into Taylor basis
dataTaylor{1} = projectDataDisc2DataTaylor(dataDisc(:,:,1), globM, globMDiscTaylor);
% Limit in Taylor basis
dataTaylor{1} = applySlopeLimiterTaylor(g, dataTaylor{1}, markV0TbdrD{1}, dataV0T{1}, basesOnQuad, type);
% Project back to original basis
dataDisc(:,:,1) = projectDataTaylor2DataDisc(dataTaylor{1}, globM, globMDiscTaylor);

% Compute Taylor cofficients of water height
dataTaylor{1} = projectDataDisc2DataTaylor(dataDisc(:,:,1) - zbDisc, globM, globMDiscTaylor);

% Project into Taylor basis
dataTaylor{2} = projectDataDisc2DataTaylor(dataDisc(:,:,2), globM, globMDiscTaylor);
dataTaylor{3} = projectDataDisc2DataTaylor(dataDisc(:,:,3), globM, globMDiscTaylor);

valCentroid{1} = dataTaylor{2}(:,1) ./ dataTaylor{1}(:,1);
valCentroid{2} = dataTaylor{3}(:,1) ./ dataTaylor{1}(:,1);

grad = { [ valCentroid{1}, dataTaylor{2}(:,2)-valCentroid{1} .* dataTaylor{1}(:,2), dataTaylor{2}(:,3)-valCentroid{1} .* dataTaylor{1}(:,3) ]; ...
         [ valCentroid{2}, dataTaylor{3}(:,2)-valCentroid{2} .* dataTaylor{1}(:,2), dataTaylor{3}(:,3)-valCentroid{2} .* dataTaylor{1}(:,3) ] };

hV0T = computeFuncDiscAtPoints(dataTaylor{1}, basesOnQuad.phiTaylorV0T(:,:,1:3));
valV0T = { computeFuncDiscAtPoints(dataTaylor{2}, basesOnQuad.phiTaylorV0T(:,:,1:3)) ./ hV0T; computeFuncDiscAtPoints(dataTaylor{3}, basesOnQuad.phiTaylorV0T(:,:,1:3)) ./ hV0T };

% Limit in Taylor basis
[alphaE, Q] = computeVectorVertexBasedLimiter(g, grad, valV0T, markV0TbdrD{2}, dataV0T(2:3), variant, angle);

rotCoef1 = Q{1,1}.^2;
rotCoef2 = Q{1,1}.*Q{1,2};
rotCoef3 = Q{2,1}.^2;
rotCoef4 = Q{2,1}.*Q{2,2};
rotCoef5 = Q{1,2}.^2;
rotCoef6 = Q{2,2}.^2;

prod1 = repmat(valCentroid{1}, [1 2]) .* dataTaylor{1}(:,2:3);
prod2 = repmat(valCentroid{2}, [1 2]) .* dataTaylor{1}(:,2:3);

diff1 = dataTaylor{2}(:,2:3) - prod1;
diff2 = dataTaylor{3}(:,2:3) - prod2;

for i = 1:2
  dataTaylor{2}(:,i+1) = prod1(:,i) + alphaE{1} .* (rotCoef1 .* diff1(:,i) + rotCoef2 .* diff2(:,i)) + alphaE{2} .* (rotCoef3 .* diff1(:,i) + rotCoef4 .* diff2(:,i));
  dataTaylor{3}(:,i+1) = prod2(:,i) + alphaE{1} .* (rotCoef2 .* diff1(:,i) + rotCoef5 .* diff2(:,i)) + alphaE{2} .* (rotCoef4 .* diff1(:,i) + rotCoef6 .* diff2(:,i));
end % for

% Project back to original basis
dataDisc(:,:,2) = projectDataTaylor2DataDisc(dataTaylor{2}, globM, globMDiscTaylor);
dataDisc(:,:,3) = projectDataTaylor2DataDisc(dataTaylor{3}, globM, globMDiscTaylor);

end % function
