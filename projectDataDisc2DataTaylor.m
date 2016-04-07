% Convert the DG/modal into a Taylor basis representation.
%
%===============================================================================
%> @file projectDataDisc2DataTaylor.m
%>
%> @brief Convert the DG/modal into a Taylor basis representation.
%===============================================================================
%>
%> @brief Convert the DG/modal into a Taylor basis representation.
%>
%> It performs an L2-projection from the DG/modal into a Taylor basis
%> representation. On a triangle @f$T_k@f$ it is defined as
%> @f[
%>   \sum_{j=1}^N C^\mathrm{DG}_{kj} 
%>               \int_{T_k} \varphi_{ki} \varphi_{kj} \mathrm{d}\mathbf{x}
%>   = \sum_{l=1}^{N_\mathrm{Taylor}} C^\mathrm{Taylor}_{kl} 
%>               \int_{T_k} \varphi_{ki} \phi_{kl} \mathrm{d}\mathbf{x} \,,
%>   \forall i \in \{1,\ldots,N\}\,.
%> @f]
%> In matrix form it can be rewritten as
%> @f[
%>   \mathsf{M} \mathbf{C}^\mathrm{DG} = 
%>      \mathsf{M}^\mathrm{DG,Taylor} \mathbf{C}^\mathrm{Taylor} \,.
%> @f]
%>
%> @param  dataDisc        Coefficient matrix @f$\mathbf{C}^\mathrm{DG}@f$ of the
%>                         DG/modal basis representation. @f$[K \times N]@f$
%> @param  globMDisc       Mass matrix @f$\mathsf{M}@f$ in the DG/modal basis.
%>                         @f$[KN \times KN]@f$
%> @param  globMDiscTaylor Projection matrix @f$\mathsf{M}^\mathrm{DG,Taylor}@f$.
%>                         @f$[KN \times KN]@f$
%>                   
%> @retval The representation matrix of the Taylor basis representation
%>         @f$\mathbf{C}^\mathrm{Taylor}@f$. @f$[K \times N]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function dataTaylor = projectDataDisc2DataTaylor(dataDisc, globMDisc, globMDiscTaylor)
[K, N] = size(dataDisc);
validateattributes(globMDisc, {'numeric'}, {'size', [K*N K*N]}, mfilename, 'globMDisc');
validateattributes(globMDiscTaylor, {'numeric'}, {'size', [K*N K*N]}, mfilename, 'globMDiscTaylor');
dataTaylor = reshape(globMDiscTaylor \ ( globMDisc * reshape(dataDisc', [K*N 1]) ), [N K])';
end % function
