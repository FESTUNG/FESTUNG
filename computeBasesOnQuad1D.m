% Evaluate tensor product basis functions and their gradients in quadrature
% points on the reference interval and stores them in a struct.

%===============================================================================
%> @file computeBasesOnQuad1D.m
%>
%> @brief Evaluate tensor product basis functions and their gradients in 
%>        quadrature points on the reference interval and stores them in a struct.
%===============================================================================
%>
%> @brief Evaluate tensor product basis functions and their gradients in 
%>        quadrature points on the reference interval and stores them in a struct.
%>
%> It evaluates the basis functions provided by <code>phi1D()</code>
%> and <code>gradPhi1D()</code> in all quadrature points for all 
%> required orders on the reference interval @f$\hat{T} = [0, 1]@f$.
%> 
%> The quadrature points on the reference interval @f$q_r@f$ are provided by 
%> <code>quadRule1D()</code>.
%>
%> All struct variables are @f$\#\mathcal{P} \times 1@f$ <code>cell</code>-arrays, 
%> with @f$\mathcal{P} = \{2p, 2p+1\}@f$ the set of required polynomial orders. 
%> Note, for @f$p=0@f$ only order @f$1@f$ is provided.
%>
%> This function computes the following struct variables (dimensions given for
%> each order):
%> - <code>phi0D</code>: @f$\hat{\phi}_i (0)@f$ and @f$\hat{\phi}_i (1)@f$ 
%>                                @f$[N \times 2]@f$
%> - <code>phi1D</code>: @f$\hat{\phi}_i (q_r) \; [R \times N]@f$
%> - <code>gradPhi1D</code>: @f$\partial_{\hat{x}}\hat{\phi}_i(q_r)\;[R \times N]@f$
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
function basesOnQuad = computeBasesOnQuad1D(p, basesOnQuad, requiredOrders)
if nargin < 3
  if p > 0
    requiredOrders = [2*p, 2*p+1]; 
  else
    requiredOrders = 1; 
  end % if
end % if
%
N = p+1;

basesOnQuad.phi0D = cell(max(requiredOrders), 1);
basesOnQuad.phi1D = cell(max(requiredOrders), 1);
basesOnQuad.gradPhi1D = cell(max(requiredOrders), 1);

for qOrd = requiredOrders
  [Q, ~] = quadRule1D(qOrd);  R = length(Q); 
  basesOnQuad.phi0D{qOrd} = zeros(N, 2);
  basesOnQuad.phi1D{qOrd} = zeros(R, N);
  basesOnQuad.gradPhi1D{qOrd} = zeros(R, N);
  for i = 1 : N
    basesOnQuad.phi1D{qOrd}(:, i) = phi1D(i, Q);
    basesOnQuad.gradPhi1D{qOrd}(:, i) = gradPhi1D(i, Q);
    for n = 1 : 2
      basesOnQuad.phi0D{qOrd}(i, n) = phi1D(i, n-1);
    end % for n
  end % for i
end % for qOrd
end  % function