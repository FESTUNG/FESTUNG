% Evaluates the i-th Taylor basis function in physical domain.
%
%===============================================================================
%> @file phiTaylorPhy.m
%>
%> @brief Evaluates the i-th Taylor basis function in physical domain.
%===============================================================================
%>
%> @brief Evaluates the @f$i@f$-th Taylor basis function in physical domain.
%>
%> The local Taylor basis on element @f$T_k@f$ is defined by
%> @f{eqnarray*}{
%> \phi_{k1} =& 1 \;, \\
%> \phi_{k2} =& \frac{x_k^1 - x_{k\mathrm{c}}^1}{\Delta (x_k^1)} \;, \\
%> \phi_{k3} =& \frac{x_k^2 - x_{k\mathrm{c}}^2}{\Delta (x_k^2)} \;, \\
%> \phi_{ki} =& \frac{(\mathbf{x} - \mathbf{x}_{k\mathrm{c}})^{\mathbf{a}_i} 
%>  - \overline{(\mathbf{x} - \mathbf{x}_{k\mathrm{c}})^{\mathbf{a}_i}}}{\mathbf{a}_i! \, (\Delta \mathbf{x}_k)^{\mathbf{a}_i}} 
%>  \quad\text{for}\quad i \ge 4 \,,
%> @f}
%> where @f$\mathbf{x}_{kc}@f$ is the centroid of triangle @f$T_k@f$, 
%> @f$\mathbf{a} = [a^1,a^2]^T@f$ is a two-dimensional multi-index,
%> and for @f$v:T_k\rightarrow \mathbb{R}@f$, let 
%> @f$\overline{v} := \frac{1}{|T_k|}\int_{T_k} v(\mathbf{x})\,\mathrm{d}\mathbf{x}@f$.
%>
%> Here, we used some standard notation for multi-indices:
%> @f{eqnarray*}{
%>  \mathbf{a} \pm \mathbf{b}  =& [a^1 \pm b^1, a^2 \pm b^2]^T \;, \\
%>  |\mathbf{a}|              :=& a^1 + a^2                   \;, \\
%>  \mathbf{a}!               :=& a^1! a^2!                   \;, \\
%>  \mathbf{x}^\mathbf{a}     :=& (x^1)^{a^1} (x^2)^{a^2}     \;, \\
%>  \partial^\mathbf{a}       :=& \partial^{|\mathbf{a}|} \big/ \partial(x^1)^{a^1}\,\partial(x^2)^{a^2} \;.
%> @f}
%> 
%> @param  g   The lists describing the geometric and topological 
%>             properties of a triangulation (see 
%>             <code>generateGridData()</code>) 
%>             @f$[1 \times 1 \text{ struct}]@f$
%> @param  i   The index of the basis function.
%> @param  X1  A list of @f$x^1@f$ coordinates.
%>             @f$[K \times n_\mathrm{Points}@f$]
%> @param  X2  A list of @f$x^2@f$ coordinates.
%>             @f$[K \times n_\mathrm{Points}@f$]
%> @retval ret The @f$i@f$-th basis function in all points specified by
%>             <code>X1</code>, <code>X2</code>. It holds <code>size(X1) ==
%>             size(X2) == size(ret)</code>.
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
function ret = phiTaylorPhy(g, i, X1, X2)
assert(size(X1, 1) == g.numT, 'At least one point required in each element')
assert(isequal(size(X1), size(X2)), 'X1 and X2 must have the same size')

% Determine quadrature rule and mapping to physical elements
qOrd = ceil((sqrt(8*i+1)-3)/2);
[Q1, Q2, W] = quadRule2D(qOrd);
Q2X1 = @(X1, X2) g.B(:, 1, 1) * X1 + g.B(:, 1, 2) * X2 + g.coordV0T(:, 1, 1) * ones(size(X1));
Q2X2 = @(X1, X2) g.B(:, 2, 1) * X1 + g.B(:, 2, 2) * X2 + g.coordV0T(:, 1, 2) * ones(size(X1));

% Extract dimensions and compute scaling parameters
R = length(W); K = g.numT; numP = size(X1, 2);
dX1 = repmat(2 ./ (max(g.coordV0T(:,:,1),[],2) - min(g.coordV0T(:,:,1),[],2)), [1 numP]);
dX2 = repmat(2 ./ (max(g.coordV0T(:,:,2),[],2) - min(g.coordV0T(:,:,2),[],2)), [1 numP]);

% Evaluate basis functions
switch i
  case 1 % (0,0)
    ret = ones(K, numP);
  case 2 % (1,0)
    ret = (X1 - repmat(g.baryT(:, 1), [1 numP])) .* dX1;
  case 3 % (0,1)
    ret = (X2 - repmat(g.baryT(:, 2), [1 numP])) .* dX2;
  case 4 % (2,0)
    ret = ( 0.5 * ( X1 - repmat(g.baryT(:, 1), [1 numP]) ).^2 - ...
          repmat( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ).^2 * W', [1 numP]) ) .* (dX1 .* dX1);
  case 5 % (1,1)
    ret = ( ( X1 - repmat(g.baryT(:, 1), [1 numP]) ) .* ( X2 - repmat(g.baryT(:, 2), [1 numP]) ) - ...
                   repmat( 2 * ( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ) .* ...
                                 ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ) ) * W', [1 numP]) ) .* (dX1 .* dX2);
  case 6 % (0,2)
    ret = ( 0.5 * ( X2 - repmat(g.baryT(:, 2), [1 numP]) ).^2 - ...
          repmat( ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ).^2 * W', [1 numP]) ) .* (dX2 .* dX2);
  case 7 % (3,0)
    ret = ( ( X1 - repmat(g.baryT(:, 1), [1 numP]) ).^3 / 6 - ...
          repmat( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ).^3 * W' / 3, [1 numP]) ) .* dX1.^3;
  case 8 % (2,1)
    ret = ( 0.5 * ( X1 - repmat(g.baryT(:, 1), [1 numP]) ).^2 .* ( X2 - repmat(g.baryT(:, 2), [1 numP]) ) - ...
          repmat( ( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ).^2 .* ...
                    ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ) ) * W', [1 numP]) ) .* (dX1.^2 .* dX2);
  case 9 % (1,2)
    ret = ( 0.5 * ( X1 - repmat(g.baryT(:, 1), [1 numP]) ) .* ( X2 - repmat(g.baryT(:, 2), [1 numP]) ).^2 - ...
          repmat( ( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ) .* ...
                    ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ).^2 ) * W', [1 numP]) ) .* (dX1 .* dX2.^2);
  case 10 % (0,3)
    ret = ( ( X2 - repmat(g.baryT(:, 2), [1 numP]) ).^3 / 6 - ...
          repmat( ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ).^3 * W' / 3, [1 numP]) ) .* dX2.^3;
  case 11 % (4,0)
    ret = ( ( X1 - repmat(g.baryT(:, 1), [1 numP]) ).^4 / 24 - ...
          repmat( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ).^4 * W' / 12, [1 numP]) ) .* dX1.^4;
  case 12 % (3,1)
    ret = ( ( X1 - repmat(g.baryT(:, 1), [1 numP]) ).^3 .* ( X2 - repmat(g.baryT(:, 2), [1 numP]) ) / 6 - ...
          repmat( ( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ).^3 .* ...
                    ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ) ) * W' / 3, [1 numP]) ) .* (dX1.^3 .* dX2);
  case 13 % (2,2)
    ret = ( ( X1 - repmat(g.baryT(:, 1), [1 numP]) ).^2 .* ( X2 - repmat(g.baryT(:, 2), [1 numP]) ).^2 / 4 - ...
          repmat( ( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ).^2 .* ...
                    ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ).^2 ) * W' / 2, [1 numP]) ) .* (dX1 .* dX2).^2;
  case 14 % (1,3)
    ret = ( ( X1 - repmat(g.baryT(:, 1), [1 numP]) ) .* ( X2 - repmat(g.baryT(:, 2), [1 numP]) ).^3 / 6 - ...
          repmat( ( ( Q2X1(Q1, Q2) - repmat(g.baryT(:, 1), [1 R]) ) .* ...
                    ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ).^3 ) * W' / 3, [1 numP]) ) .* (dX1 .* dX2.^3);
  case 15 % (0,4)
    ret = ( ( X2 - repmat(g.baryT(:, 2), [1 numP]) ).^4 / 24 - ...
          repmat( ( Q2X2(Q1, Q2) - repmat(g.baryT(:, 2), [1 R]) ).^4 * W' / 12, [1 numP]) ) .* dX2.^4;
  otherwise
    error('Taylor basis functions for p>4 not implemented')
end % switch
end

