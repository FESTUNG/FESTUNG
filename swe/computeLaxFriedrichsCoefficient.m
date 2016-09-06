% Computes the Lax-Friedrichs flux coefficient for the 2D Shallow-Water 
% Equations.

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
%> @f[ \begin{bmatrix}
%> 0 & \nu_x & \nu_y \\
%> gh\nu_x - u^2\nu_x - uv\nu_y & 2u\nu_x + v\nu_y & u\nu_y \\
%> gh\nu_y - uv\nu_x - v^2\nu_y & v\nu_x & u\nu_x + 2v\nu_y
%> \end{bmatrix} @f].
%>
%> Using a formal computation toolbox, one can see that the eigenvalues are
%> 
%> @f$ \lambda_1 = u_n - a,~\lambda_2 = u_n,~\lambda_3 = u_n + a,@f$
%> with @f$ u_n = u \nu_x + v \nu_y,@f$ and @f$a = \sqrt{gh}@f$, where @f$a > 0@f$ 
%> must hold.
%> The largest absolute value @f$ \lambda = \max( |{\lambda_1}|,|{\lambda_2}|,|{\lambda_3}| ) @f$ 
%> equals @f$ |{u_n}| + a@f$, as can be seen by the following distinction of 
%> cases:
%>
%> if @f$ 0 \leq u_n \leq a@f$, then @f$ \lambda = \max( -u_n+a, u_n, u_n + a ) = |{u_n}| + a,@f$ \n
%> if @f$ 0 < a < u_n@f$, then @f$ \lambda = \max( u_n-a, u_n, u_n + a ) = |{u_n}| + a,@f$ \n
%> if @f$ u_n < 0 < |{u_n}| \leq a@f$, then @f$ \lambda = \max( -u_n+a, u_n, u_n + a ) = |{u_n}| + a,@f$ \n
%> if @f$ u_n < 0 < a < |{u_n}|@f$, then @f$ \lambda = \max( -u_n+a, u_n, -u_n - a ) = |{u_n}| + a,@f$
%>
%> @param cQ0E0T    The averaged values of the height and velocities as computed
%>                  by e.g. computeAveragedVariablesQ0E0Tint() in each 
%>                  quadrature point of each local edge of a particular local 
%>                  index. @f$[3 \times 1 \text{ cell}]@f$
%> @param nuQ0E0T   The normal components for each quadrature point of each 
%>                  local edge of a particular local index. @f$[2 \times 1 \text{ cell}]@f$
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

validateattributes(cQ0E0T, {'cell'}, {'size', [3 1]}, mfilename, 'cQ0E0T');
validateattributes(nuQ0E0T, {'cell'}, {'size', [1 2]}, mfilename, 'nuQ0E0T');

assert(isscalar(gConst) && gConst >= 0, 'Gravitational constant has to be a non-negative scalar.')

lambda = setNaN2Zero(abs(cQ0E0T{2} .* nuQ0E0T{1} + cQ0E0T{3} .* nuQ0E0T{2})) + sqrt(gConst * cQ0E0T{1});
end % function
