% Applies a sequential slope limiter to a discrete function and products of
% functions with the former, all given in (modal) DG basis.

%===============================================================================
%> @file
%>
%> @brief Applies a sequential slope limiter to a discrete function and 
%>        products of functions with the former, all given in (modal) DG
%>        basis.
%===============================================================================
%>
%> @brief Applies a sequential slope limiter to a discrete function and 
%>        products of functions with the former, all given in (modal) DG
%>        basis.
%>
%> It applies the sequential slope limiting operator 
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
%> (cf. <code>computeSequentialVertexBasedLimiter()</code>).
%>
%> @param  g              The lists describing the geometric and topological 
%>                        properties of a triangulation (see 
%>                        <code>generateGridData()</code>) 
%>                        @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc       The representation matrix of the unlimited functions 
%>                        @f$H_h@f$ and @f$\mathbf{u}_h@f$.
%> @param  zbDisc         The representation matrix of the problem
%>                        bathymetry @f$b_h@f$ @f$[K \times N]@f$.
%> @param  markV0TbdrD    <code>cell</code> that for each unknown marks
%>                        each triangles (Dirichlet boundary) vertices on 
%>                        which additional function values should be
%>                        considered during the  slope limiting routine. 
%>                        @f$[3 \times 1]@f$
%> @param  dataV0T        The function values for each unknown for 
%>                        (Dirichlet boundary) vertices specified by 
%>                        <code>markV0TbdrD</code>. @f$[3 \times 1]@f$
%> @param  globM          The mass matrix @f$\mathsf{M}@f$ in DG basis.
%>                        @f$[KN \times KN]@f$.
%> @param globMDiscTaylor The transformation matrix @f$\mathsf{M}^\mathrm{DG,Taylor}@f$.
%>                        @f$[KN \times KN]@f$.
%> @param  basesOnQuad    A struct containing precomputed values of (Taylor) basis
%>                        functions on quadrature points. Must provide at
%>                        least phiTaylorV0T.
%> @param  type           The type of slope limiter to be used. [<code>string</code>]
%> @retval dataDisc       The representation matrix of the limited function
%>                        @f$\mathsf{\Phi}c_h@f$. @f$[K \times N \times d+1]@f$
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
function dataDisc = applySequentialSlopeLimiterDisc(g, dataDisc, zbDisc, markV0TbdrD, dataV0T, globM, globMDiscTaylor, basesOnQuad, type)

sizeUnknowns = size(dataDisc);
dataTaylor = zeros(sizeUnknowns);
% Project into Taylor basis
dataTaylor(:,:,1) = projectDataDisc2DataTaylor(dataDisc(:,:,1), globM, globMDiscTaylor);
% Limit in Taylor basis
dataTaylor(:,:,1) = applySlopeLimiterTaylor(g, dataTaylor(:,:,1), markV0TbdrD{1}, dataV0T{1}, basesOnQuad, type);
% Project back to original basis
dataDisc(:,:,1) = projectDataTaylor2DataDisc(dataTaylor(:,:,1), globM, globMDiscTaylor);

% Compute Taylor cofficients of water height
dataTaylor(:,:,1) = projectDataDisc2DataTaylor(dataDisc(:,:,1) - zbDisc, globM, globMDiscTaylor);

hV0T = computeFuncDiscAtPoints(dataTaylor(:,1:3,1), basesOnQuad.phiTaylorV0T(:,:,1:3));

for i = 2 : sizeUnknowns(3)
  % Project into Taylor basis
  dataTaylor(:,:,i) = projectDataDisc2DataTaylor(dataDisc(:,:,i), globM, globMDiscTaylor);
  
  % prepare sequential limiting
  valCentroid = dataTaylor(:,1,i) ./ dataTaylor(:,1,1);
  prod = repmat(valCentroid, [1 2]) .* dataTaylor(:,2:3,1);
  % compute Limiting factor
  alphaE = computeVertexBasedLimiter(g, valCentroid, computeFuncDiscAtPoints(dataTaylor(:,1:3,i), basesOnQuad.phiTaylorV0T(:,:,1:3)) ./ hV0T, markV0TbdrD{i}, dataV0T{i}./dataV0T{1});
  
  % Limit in Taylor basis
  dataTaylor(:,2:3,i) = prod + repmat(alphaE, [1 2]) .* (dataTaylor(:,2:3,i) - prod);
  % Project back to original basis
  dataDisc(:,:,i) = projectDataTaylor2DataDisc(dataTaylor(:,:,i), globM, globMDiscTaylor);
end % for

end % function
