% Provides weights and time levels for explicit strong stability preserving (SSP)
% Runge-Kutta methods up to order 3.

%===============================================================================
%> @file
%>
%> @brief Provides weights and time levels for explicit strong stability 
%>        preserving (SSP) Runge-Kutta methods up to order 3.
%===============================================================================
%>
%> @brief Provides weights and time levels for explicit strong stability 
%>        preserving (SSP) Runge-Kutta methods up to order 3.
%>
%> Given a time-dependent system 
%> @f[
%> \partial_t \mathbf{C}(t) = \mathbf{S}(\mathbf{C}(t),t), \quad t \in J = (0,t^\mathrm{end}).
%> @f]
%> Let @f$0 = t^1 < t^2 < \ldots < t^\mathrm{end}@f$ be a not necessarily
%> equidistant decomposition of the time interval @f$J@f$ and let 
%> @f$\Delta t^n = t^{n+1} - t^{n}@f$ denote the time step size.
%> The update-scheme of the @f$s@f$-step Runge-Kutta method (of order @f$s@f$)
%> is then defined as
%> @f[
%>  \mathbf{C}^{(0)} := \mathbf{C}^n,
%> @f]@f[
%>  \mathbf{C}^{(i)} := \omega_i \mathbf{C}^n + (1-\omega_i) (\mathbf{C}^{(i-1)} + \Delta t^n \mathbf{S}(\mathbf{C}^{(i-1)}, t^{(i)}), \quad \mathrm{for} \quad i = 1,\ldots,s
%> @f]@f[
%>  \mathbf{C}^{n+1} := \mathbf{C}^{(s)}
%> @f]
%> with @f$\mathbf{C}^n = \mathbf{C}(t^n)@f$.
%>
%> The Runge-Kutta methods are total variation diminishing (TVD).
%> For details refer to
%> Gottlieb, Sigal, and Chi-Wang Shu. "Total Variation Diminishing 
%> Runge-Kutta Schemes". Mathematics of Computation 67.221 (1998): 73–85.
%> doi: 10.1090/S0025-5718-98-00913-2.
%>
%> @param  ord   The order of the time stepping scheme.
%> @param  tau   The time step size.
%> @param  t0    The current time-level @f$t^n@f$.
%> @retval t     The vector of time-levels @f$t^{(i)}@f$. @f$[1 \times s]@f$
%> @retval omega The associated weights. @f$[1 \times s]@f$
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
function [t, omega] = rungeKuttaExplicit(ord, tau, t0)
switch ord
  case 1
    omega = 0;
    t = t0;
  case 2
    omega = [0, 0.5];
    t = t0 + [0, 1] * tau;
  case 3
    omega = [0, 3/4, 1/3];
    t = t0 + [0, 1, 0.5] * tau;
  otherwise
    error('Order %d not implemented', ord);
end
end % function