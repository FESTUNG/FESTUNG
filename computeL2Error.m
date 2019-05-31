% Compute the L2-error of a DG/modal basis representation for a given 
% analytical solution.

%===============================================================================
%> @file
%>
%> @brief Compute the L2-error of a DG/modal basis representation for a given 
%>        analytical solution.
%===============================================================================
%>
%> @brief Compute the L2-error of a DG/modal basis representation for a given 
%>        analytical solution.
%>
%> The discretization error @f$\|c_h(t)-c(t)\|_{L^2(\Omega)}@f$ at time @f$t@f$
%> gives the @f$L^2@f$-norm of the difference between the discrete solution 
%> @f$c_h(t)@f$ and the analytical solution @f$c(t)@f$.
%>
%> To compute the discretization error, we must evaluate
%> @f[
%>  \|c_h(t)-c(t)\|^2_{L^2(\Omega)} = 
%>    \sum_{T_k\in\mathcal{T}_h} \int_{T_k} (c_h(t)-c(t))^2 \mathrm{d}\mathbf{x}.
%> @f]
%> The DG/modal basis representation of the discrete solution @f$c(t)@f$,
%> is given as
%> @f[
%>   c_h(t, \mathbf{x}) = \sum_{j=1}^N C_{kj}(t) \varphi_{kj}(\mathbf{x}) \,.
%> @f]
%> Using an affine mapping 
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$ from 
%> the reference triangle @f$\hat{T} = \{(0,0), (1,0), (0,1)\}@f$, defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,,
%> @f]
%> we obtain
%> @f[
%>   \|c_h(t)-c(t)\|^2_{L^2(\Omega)} = 
%>    2 \sum_{T_k\in\mathcal{T}_h}  |T_k| \int_{\hat{T}}
%>    (\sum_{j=1}^N C_{kj}(t) \hat{\varphi}_j(\hat{\mathbf{x}}) - c \circ
%>    \mathbf{F}_k(\hat{\mathbf{x}}))^2 \mathrm{d} \hat{\mathbf{x}} \,.
%> @f]
%> Applying a quadrature rule as given by <code>quadRule2D()</code>, we
%> obtain
%> @f[
%>   \|c_h(t)-c(t)\|^2_{L^2(\Omega)} \approx
%>    2 \sum_{T_k\in\mathcal{T}_h} |T_k| \sum_{r=1}^R \omega_r
%>      (\sum_{j=1}^N C_{kj}(t) \hat{\varphi}_j(\hat{\mathbf{q}}_r) - c \circ
%>      \mathbf{F}_k(\hat{\mathbf{q}}_r))^2 =
%>    2 \begin{bmatrix} |T_1| \\ \vdots \\ |T_K| \end{bmatrix} \cdot \left(
%>    \begin{bmatrix} C_{11}&\dots&C_{1N} \\ \vdots&&\vdots \\ C_{K1}&\dots&C_{KN}
%>    \end{bmatrix} \begin{bmatrix}
%>    \hat{\varphi}_1(\hat{\mathbf{q}}_1)&\dots&\hat{\varphi}_1(\hat{\mathbf{q}}_R)\\
%>    \vdots & & \vdots \\
%>    \hat{\varphi}_N(\hat{\mathbf{q}}_1)&\dots&\hat{\varphi}_N(\hat{\mathbf{q}}_R)
%>    \end{bmatrix}-c(\mathsf{X}^1, \mathsf{X}^2) \right)^2
%>    \begin{bmatrix} \omega_1 \\ \vdots \\ \omega_R \end{bmatrix} \,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc   The representation matrix of the DG/modal basis
%>                    representation @f$[K \times N]@f$
%> @param  funcCont   A function handle for the continuous solution
%> @param  qOrd       The order of the quadrature rule provided by 
%>                    <code>quadRule2D()</code>
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret The discretization error @f$[\text{scalar}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function err = computeL2Error(g, dataDisc, funcCont, qOrd, basesOnQuad)
% Check function arguments that are directly used
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(dataDisc, {'numeric'}, {'size', [g.numT NaN]}, mfilename, 'dataDisc');

% Determine quadrature rule and physical coordinates
qOrd = max(qOrd,1); [Q1, Q2, W] = quadRule2D(qOrd);
R = length(W); N = size(dataDisc,2);
X1 = kron(g.B(:,1,1), Q1) + kron(g.B(:,1,2), Q2) + kron(g.coordV0T(:,1,1), ones(1,R));
X2 = kron(g.B(:,2,1), Q1) + kron(g.B(:,2,2), Q2) + kron(g.coordV0T(:,1,2), ones(1,R));

% Evaluate analytical and discrete function
cExOnQuadPts = funcCont(X1, X2); % [K x R]
cApprxOnQuadPts = dataDisc * basesOnQuad.phi2D{qOrd}(:, 1:N).'; % [K x R] = [K x N] * [N x R]

% Compute error
err = sqrt(2 * dot((cApprxOnQuadPts - cExOnQuadPts).^2 * W.', g.areaT)); 
end % function
