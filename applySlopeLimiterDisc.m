% Applies a slope limiter to a discrete function given in (modal) DG basis.

%===============================================================================
%> @file applySlopeLimiterDisc.m
%>
%> @brief Applies a slope limiter to a discrete function given in (modal) DG basis.
%===============================================================================
%>
%> @brief Applies a slope limiter to a discrete function given in (modal) DG basis.
%>
%> It applies the slope limiting operator 
%> @f[
%>  \mathsf{\Phi} := \mathsf{M}^{-1} \mathsf{M}^\mathrm{DG,Taylor} 
%>  \mathsf{\Phi}^\mathrm{Taylor} \left(\mathsf{M}^\mathrm{DG,Taylor}\right)^{-1} \mathsf{M}
%> @f]
%> to a discrete function @f$c_h: \Omega \rightarrow \mathbb{R}@f$ given as
%> representation matrix, e.g., obtained by L2-projection of a continuous
%> function (cf. <code>projectFunctCont2DataDisc()</code>).
%> @f$\mathsf{M}@f$ is the mass matrix in DG basis (cf.
%> <code>assembleMatElemPhiPhi()</code>),
%> @f$\mathsf{M}^\mathrm{DG,Taylor}@f$ is the transformation matrix from
%> (modal) DG to Taylor basis (cf.
%> <code>assembleMatElemPhiDGPhiTaylor()</code>), and
%> @f$\mathsf{\Phi}^\mathrm{Taylor}@f$ is the slope limiting operator in
%> Taylor basis (cf. <code>applySlopeLimiterTaylor()</code>).
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc   The representation matrix of the unlimited function 
%>                    @f$c_h@f$. @f$[K \times N]@f$
%> @param  markV0TbdrD <code>logical</code> arrays that mark each triangles
%>                    (Dirichlet boundary) vertices on which additional
%>                    function values should be considered during the 
%>                    slope limiting routine. @f$[K \times 3]@f$
%> @param  dataV0T    The function values for (Dirichlet boundary) vertices
%>                    specified by <code>markV0TbdrD</code>. @f$[K \times 3]@f$
%> @param  globM      The mass matrix @f$\mathsf{M}@f$ in DG basis.
%>                    @f$[KN \times KN]@f$.
%> @param globMDiscTaylor  The transformation matrix @f$\mathsf{M}^\mathrm{DG,Taylor}@f$.
%>                         @f$[KN \times KN]@f$.
%> @param  basesOnQuad  A struct containing precomputed values of (Taylor) basis
%>                      functions on quadrature points. Must provide at
%>                      least phiTaylorV0T.
%> @param  type       The type of slope limiter to be used. [<code>string</code>]
%> @retval dataDisc   The representation matrix of the limited function
%>                    @f$\mathsf{\Phi}c_h@f$. @f$[K \times N]@f$
%> @retval  minMaxV0T  Two matrices with minimum or maximum centroid values,
%>                     respectively, of the patch of elements surrounding each
%>                     vertex of each element as computed by 
%>                     <code>computeMinMaxV0TElementPatch()</code>
%>                     @f$[2 \times 1 \mathrm{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>                      Modified 09/02/16 by Hennes Hajduk
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
function [dataDisc, minMaxV0T] = applySlopeLimiterDisc(g, dataDisc, markV0TbdrD, dataV0T, globM, globMDiscTaylor, basesOnQuad, type)
% Project into Taylor basis
dataTaylor = projectDataDisc2DataTaylor(dataDisc, globM, globMDiscTaylor);

% Limit in Taylor basis
[dataTaylor, minMaxV0T] = applySlopeLimiterTaylor(g, dataTaylor, markV0TbdrD, dataV0T, basesOnQuad, type);

% Project back to original basis
dataDisc = projectDataTaylor2DataDisc(dataTaylor, globM, globMDiscTaylor);
end

