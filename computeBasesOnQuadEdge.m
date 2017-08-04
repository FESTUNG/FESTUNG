% Evaluate edge basis functions in quadrature points on the
% reference edge and stores them in a struct.

%===============================================================================
%> @file computeBasesOnQuadEdge.m
%>
%> @brief Evaluate edge basis functions in quadrature points on the
%> reference edge and stores them in a struct.
%===============================================================================
%>
%> @brief Evaluate edge basis functions in quadrature points on the
%> reference edge and stores them in a struct.
%>
%> It evaluates the edge basis functions provided by <code>phi1D()</code> 
%> in all quadrature points for all required orders 
%> on the reference edge @f$\hat{E} = [0,1]@f$.
%> 
%> The quadrature points @f$q_r@f$ on the reference edge @f$[0,1]@f$ are 
%> provided by <code>quadRule1D()</code>.
%> 
%> All struct variables are @f$\#\mathcal{P} \times 1@f$ <code>cell</code>-arrays, 
%> with @f$\mathcal{P} = \{2p, 2p+1\}@f$ the set of required polynomial orders. 
%> Note, for @f$p=0@f$ only order @f$1@f$ is provided.
%> 
%> This function computes the following struct variables (dimensions given for
%> each order):
%> - <code>mu</code>: @f$\hat{\mu}_i (q_r) \; [R \times \bar{N}]@f$
%> - <code>thetaMu</code>: @f$\hat{\mu}_i \circ \hat{\beta}_{kn}(q_r)
%>                         \; [R \times \bar{N} \times 2]@f$
%> 
%> Here, the mapping @f$\hat{\beta}_{kn}@f$ adapts the edge orientation to
%> match the definitions of edge basis function from an element view and an
%> edge view. It is defined as
%> 
%> @f[
%>  \hat{\beta}_{kn} =
%>  \begin{cases}
%>    s   &\text{if } T_k \text{ is the first element in the edge-local indexing}, \\
%>    1-s &\text{if } T_k \text{ is the second element in the edge-local indexing}\,.
%>  \end{cases}
%> @f]
%>
%> @param  N          The number of local edge degrees of freedom. For polynomial
%>                    order @f$p@f$, it is given as @f$\bar{N} = p+1@f$
%>                    @f$[\text{scalar}]@f$
%> @param  basesOnQuadEdge A (possibly empty) struct to which the computed
%>                    arrays are added. @f$[\text{struct}]@f$
%> @param  requiredOrders (optional) An array providing a list of all
%>                    required quadrature orders. @f$[\text{scalar}]@f$
%>
%> @retval  basesOnQuadEdge A struct with the computed array.  @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
function basesOnQuadEdge = computeBasesOnQuadEdge(N, basesOnQuadEdge, requiredOrders)
validateattributes(basesOnQuadEdge, {'struct'}, {}, mfilename, 'basesOnQuadEdge')

if nargin < 3
    p = N - 1;
    if p > 0
        requiredOrders = [2*p, 2*p+1];
    else
        requiredOrders = 1;
    end % if
end % if

% Precompute basis functions
basesOnQuadEdge.mu = cell(max(requiredOrders),1);
basesOnQuadEdge.thetaMu = cell(max(requiredOrders),1);
for it = 1 : length(requiredOrders)
    ord = requiredOrders(it);
    [Q, ~] = quadRule1D(ord);
    R = length(Q);
    
    basesOnQuadEdge.mu{ord} = zeros(R, N);
    for i = 1 : N
        basesOnQuadEdge.mu{ord}(:, i) = phi1D(i, Q);
    end
    
    basesOnQuadEdge.thetaMu{ord} = zeros(R, N, 2);
    basesOnQuadEdge.thetaMu{ord}(:, :, 1) = basesOnQuadEdge.mu{ord};
    basesOnQuadEdge.thetaMu{ord}(:, :, 2) = flipud(basesOnQuadEdge.mu{ord});
end % for
end % function
