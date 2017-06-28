% Evaluate tensor product basis functions and their gradients in quadrature
%  points on the reference square and stores them in a struct.

%===============================================================================
%> @file computeBasesOnQuadTensorProduct.m
%>
%> @brief Evaluate tensor product basis functions and their gradients in 
%>        quadrature points on the reference square and stores them in a struct.
%===============================================================================
%>
%> @brief Evaluate tensor product basis functions and their gradients in 
%>        quadrature points on the reference square and stores them in a struct.
%>
%> It evaluates the basis functions provided by <code>phiTensorProduct()</code>
%> and <code>gradPhiTensorProduct()</code> in all quadrature points for all 
%> required orders on the reference square 
%> @f$\hat{T} = \{(0,0), (1,0), (1,1), (0,1) \}@f$
%> 
%> The quadrature points on the reference square @f$\mathbf{q}_r@f$ and the
%> unit interval @f$[0,1]@f$, @f$q_r@f$, are provided by 
%> <code>quadRule1D()</code> and its tensor product.
%>
%> All struct variables are @f$\#\mathcal{P} \times 1@f$ <code>cell</code>-arrays, 
%> with @f$\mathcal{P} = \{2p, 2p+1\}@f$ the set of required polynomial orders. 
%> Note, for @f$p=0@f$ only order @f$1@f$ is provided.
%>
%> This function computes the following struct variables (dimensions given for
%> each order):
%> - <code>phi1D</code>: @f$\hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(q_r)
%>                           \; [R \times N \times 4]@f$
%> - <code>phi2D</code>: @f$\hat{\varphi}_i(\mathbf{q}_r) \; [R \times N]@f$
%> - <code>gradPhi2D</code>: @f$\hat{\nabla} \hat{\varphi}_i(\mathbf{q}_r)
%>                               \; [R \times N \times 2]@f$
%> 
%> @param  p           The polynomial degree.
%> @param  basesOnQuad A (possibly empty) struct to which the computed
%>                     arrays are added. @f$[\text{struct}]@f$
%> @param  requiredOrders (optional) An array providing a list of all
%>                     required quadrature orders.
%>
%> @retval  basesOnQuad A struct with the computed arrays.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function basesOnQuad = computeBasesOnQuadTensorProduct(p, basesOnQuad, requiredOrders)
if nargin < 3
  if p > 0
    requiredOrders = [2*p, 2*p+1]; 
  else
    requiredOrders = 1; 
  end % if
end % if
%
N = (p+1) * (p+1);

basesOnQuad.phi1D = cell(max(requiredOrders), 1);
basesOnQuad.phi2D = cell(max(requiredOrders), 1);
basesOnQuad.gradPhi2D = cell(max(requiredOrders), 1);

for qOrd = requiredOrders
  [Q, ~] = quadRule1D(qOrd);
  [Q1, Q2, ~] = quadRuleTensorProduct(qOrd);
  R1D = length(Q); R2D = length(Q1);

  basesOnQuad.phi1D{qOrd} = zeros(R1D, N, 4);
  basesOnQuad.phi2D{qOrd} = zeros(R2D, N);
  basesOnQuad.gradPhi2D{qOrd} = zeros(R2D, N, 2);
  for i = 1 : N
    basesOnQuad.phi2D{qOrd}(:,i) = phiTensorProduct(i, Q1, Q2, @phi1D, @phi1D);
    for m = 1 : 2
      basesOnQuad.gradPhi2D{qOrd}(:,i,m) = gradPhiTensorProduct(i, m, Q1, Q2, @phi1D, @phi1D, @gradPhi1D, @gradPhi1D);
    end % for m
    for n = 1 : 4
      [QS1, QS2] = gammaMapTetra(n, Q);
      basesOnQuad.phi1D{qOrd}(:,i,n) = phiTensorProduct(i, QS1, QS2, @phi1D, @phi1D);
    end % for n
  end  % for i
end % for qOrd
end  % function