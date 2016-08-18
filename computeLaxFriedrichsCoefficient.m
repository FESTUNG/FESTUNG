% Computes the Lax-Friedrichs flux coefficient for the 2D Shallow-Water 
% Equations.
%
%===============================================================================
%> @file computeLaxFriedrichsCoefficient.m
%>
%> @brief Computes the Lax-Friedrichs flux coefficient for the 2D Shallow-
%>        Water Equations.
%===============================================================================
%>
%> @brief Computes the Lax-Friedrichs flux coefficient for the 2D Shallow-
%>        Water Equations.
%>
%> Computes the largest in absolute value of the eigenvalues of the following 
%> matrix, which is the formal Jacobian of the advective term derived by the 
%> unknowns:
%>
%> @f[
%> 0 & n_x & n_y \\
%> gHn_x - u^2n_x - uvn_y & 2un_x + vn_y & un_y \\
%> gHn_y - uvn_x - v^2n_y & vn_x & un_x + 2vn_y
%> @f].
%>
%> Using a formal computation toolbox, one can see that the eigenvalues are
%> 
%> @f$ l1 = u_n - a,~l2 = u_n,~l3 = u_n + a,@f$
%> with @f$u_n = u*nu_x + v*nu_y,@f$ and @f$a = sqrt(g*H), where @f$a > 0@f$ 
%> must hold.
%>
%> The largest absolute value @f$ \lambda = max( abs(l1),abs(l2),abs(l3) ) @f$ 
%> equals @f$ \abs(u_n) + a@f$, as can be seen by the following distinction of 
%> cases:
%>
%> if @f$ 0 <= u_n <= a, then @f$ \lambda = max( -u_n+a, u_n, u_n + a ) = abs(u_n) + a,@f$
%> if @f$ 0 < a < u_n, then @f$ \lambda = max( u_n-a, u_n, u_n + a ) = abs(u_n) + a,@f$
%> if @f$ u_n < 0 < abs(u_n) <= a, then @f$ \lambda = max( -u_n+a, u_n, u_n + a ) = abs(u_n) + a,@f$
%> if @f$ u_n < 0 < a < abs(u_n), then @f$ \lambda = max( -u_n+a, u_n, -u_n - a ) = abs(u_n) + a,@f$
%>
%> @param cQ0E0T    The averaged values of the height and velocities as computed
%>                  by e.g. computeAveragedVariablesQ0E0Tint() in each 
%>                  quadrature point of each local edge of a particular local 
%>                  index. @f${3 \times 1}@f$
%> @param nuQ0E0T   The normal components for each quadrature point of each 
%>                  local edge of a particular local index. @f${2 \times 1}@f$
%> @param gConst    The scalar valued gravitational constant.
%>
%> @retval lambda   The largest in absolute value Eigenvalue of the Jacobian of 
%>                  the advective term.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>                      Modified 08/18/16 by Hennes Hajduk
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
%>
function lambda = computeLaxFriedrichsCoefficient(cQ0E0T, nuQ0E0T, gConst)
lambda = setNaN2Zero(abs(cQ0E0T{2} .* nuQ0E0T{1} + cQ0E0T{3} .* nuQ0E0T{2})) + sqrt(gConst * cQ0E0T{1});
end % function
