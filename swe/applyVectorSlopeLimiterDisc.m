% Applies a slope limiter to a two component vector of discrete functions 
% given in (modal) DG basis.

%===============================================================================
%> @file
%>
%> @brief Applies a slope limiter to a two component vector of discrete 
%>        functions given in (modal) DG basis.
%===============================================================================
%>
%> @brief Applies a slope limiter to a two component vector of discrete 
%>        functions given in (modal) DG basis.
%>
%> It applies the slope limiting operator for a two component vector field
%> @f[
%>  \mathsf{\Phi} := \mathsf{M}^{-1} \mathsf{M}^\mathrm{DG,Taylor} 
%>  \mathsf{\Phi}^\mathrm{Taylor} \left(\mathsf{M}^\mathrm{DG,Taylor}\right)^{-1} \mathsf{M}
%> @f]
%> to a two component vector of discrete functions 
%> @f$\mathbf{u}_h: \Omega \rightarrow \mathbb{R}^2@f$ given as
%> representation matrices, e.g., obtained by L2-projection of a continuous
%> vector field (cf. <code>projectFunctCont2DataDisc()</code>).
%> @f$\mathsf{M}@f$ is the mass matrix in DG basis (cf.
%> <code>assembleMatElemPhiPhi()</code>),
%> @f$\mathsf{M}^\mathrm{DG,Taylor}@f$ is the transformation matrix from
%> (modal) DG to Taylor basis (cf.
%> <code>assembleMatElemPhiDGPhiTaylor()</code>), and
%> @f$\mathsf{\Phi}^\mathrm{Taylor}@f$ is the slope limiting operator in
%> Taylor basis (cf. <code>applySlopeLimiterTaylor()</code>) for a two 
%> component vector field.
%>
%> @param  g                The lists describing the geometric and topological 
%>                          properties of a triangulation (see 
%>                          <code>generateGridData()</code>) 
%>                          @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc         The representation matrices of the unlimited vector 
%>                          field  @f$\mathbf{u}_h@f$. @f$[K \times N \times 2]@f$
%> @param  zbDisc           The representation matrix of the problem
%>                          bathymetry @f$b_h@f$ @f$[K \times N]@f$.
%> @param  markV0TbdrD      <code>logical</code> arrays that mark each triangles
%>                          (Dirichlet boundary) vertices on which additional
%>                          function values should be considered during the 
%>                          slope limiting routine. @f$[K \times 3]@f$
%> @param  dataV0T          The vector field values for (Dirichlet boundary) vertices
%>                          specified by <code>markV0TbdrD</code>.
%>                          @f$[2 \times 1 \mathrm{cell}]@f$
%> @param  globM            The mass matrix @f$\mathsf{M}@f$ in DG basis.
%>                          @f$[KN \times KN]@f$.
%> @param  globMDiscTaylor  The transformation matrix @f$\mathsf{M}^\mathrm{DG,Taylor}@f$.
%>                          @f$[KN \times KN]@f$.
%> @param  basesOnQuad      A struct containing precomputed values of (Taylor) basis
%>                          functions on quadrature points. Must provide at
%>                          least phiTaylorV0T.
%> @param  vectorLimType    The type of vector limiter to be used. [<code>string</code>]
%> @param  variant          The variant that is used to compute rotation
%>                          matrices for directional vector limiting.
%> @param  angle            The rotation angle to be used in case a fixed
%>                          transformation is desired to compute the limiting directions.
%> @param  dataDisc         The representation matrices of the limited vector 
%>                          field  @f$\mathbf{u}_h@f$. @f$[K \times N \times 2]@f$
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
function dataDisc = applyVectorSlopeLimiterDisc(g, dataDisc, zbDisc, markV0TbdrD, dataV0T, globM, globMDiscTaylor, basesOnQuad, vectorLimType, variant, angle)

% Project into Taylor basis
dataTaylor{1} = projectDataDisc2DataTaylor(dataDisc(:,:,2), globM, globMDiscTaylor);
dataTaylor{2} = projectDataDisc2DataTaylor(dataDisc(:,:,3), globM, globMDiscTaylor);

dataTaylorLim = cell(2,1);

% Limit in Taylor basis
switch vectorLimType
  case 'direction'
    dataTaylorLim{1}(:,1) = dataTaylor{1}(:,1);  dataTaylorLim{2}(:,1) = dataTaylor{2}(:,1);
    
    valV0T = { computeFuncDiscAtPoints(dataTaylor{1}(:,1:3), basesOnQuad.phiTaylorV0T(:,:,1:3)); computeFuncDiscAtPoints(dataTaylor{2}(:,1:3), basesOnQuad.phiTaylorV0T(:,:,1:3)) };
    
    [alphaE, Q] = computeVectorVertexBasedLimiter(g, dataTaylor, valV0T, markV0TbdrD, dataV0T, variant, angle);
    
    rotCoef1 = Q{1,1}.^2;
    rotCoef2 = Q{1,1}.*Q{1,2};
    rotCoef3 = Q{2,1}.^2;
    rotCoef4 = Q{2,1}.*Q{2,2};
    rotCoef5 = Q{1,2}.^2;
    rotCoef6 = Q{2,2}.^2;
    
    for i = 2:3
      dataTaylorLim{1}(:,i) = alphaE{1} .* ( rotCoef1.*dataTaylor{1}(:,i) + rotCoef2.*dataTaylor{2}(:,i) ) + alphaE{2} .* ( rotCoef3.*dataTaylor{1}(:,i) + rotCoef4.*dataTaylor{2}(:,i) );
      dataTaylorLim{2}(:,i) = alphaE{1} .* ( rotCoef2.*dataTaylor{1}(:,i) + rotCoef5.*dataTaylor{2}(:,i) ) + alphaE{2} .* ( rotCoef4.*dataTaylor{1}(:,i) + rotCoef6.*dataTaylor{2}(:,i) );
    end % for
    
  case 'norm-momentum'
    valV0T = { computeFuncDiscAtPoints(dataTaylor{1}(:, 1:3), basesOnQuad.phiTaylorV0T(:,:,1:3)); computeFuncDiscAtPoints(dataTaylor{2}(:, 1:3), basesOnQuad.phiTaylorV0T(:,:,1:3)) };
    dataTaylorLim = applyLengthVertexBasedLimiter(g, dataTaylor, valV0T, markV0TbdrD, dataV0T);
  case 'norm-velocity'
    hV0T = projectDataDisc2DataLagr(dataDisc(:,:,1) - zbDisc);
    valV0T = { projectDataDisc2DataLagr(dataDisc(:,:,2)) ./ hV0T; projectDataDisc2DataLagr(dataDisc(:,:,3)) ./ hV0T };
    dataTaylorLim = applyLengthVertexBasedLimiter(g, dataTaylor, valV0T, markV0TbdrD, dataV0T);
  otherwise
    error('Unsupported type of vector limiting.')
end % switch

% Project back to original basis
dataDisc(:,:,2) = projectDataTaylor2DataDisc(dataTaylorLim{1}, globM, globMDiscTaylor);
dataDisc(:,:,3) = projectDataTaylor2DataDisc(dataTaylorLim{2}, globM, globMDiscTaylor);

end % function
