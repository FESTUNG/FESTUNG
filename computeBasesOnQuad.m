% Evaluate basis functions and their gradients in quadrature points on the
% reference triangle and store them in global variables.
%
%===============================================================================
%> @file computeBasesOnQuad.m
%>
%> @brief Evaluate basis functions and their gradients in quadrature points on
%>        the reference triangle and store them in global variables.
%===============================================================================
%>
%> @brief Evaluate basis functions and their gradients in quadrature points on
%>        reference triangle and store them in global variables.
%>
%> It evaluates the basis functions provided by <code>phi()</code> and
%> <code>gradPhi()</code> in all quadrature points for all required orders 
%> on the reference triangle @f$\hat{T} = \{(0,0), (1,0), (0,1) \}@f$
%> 
%> The quadrature points on the reference triangle @f$\mathbf{q}_r@f$ and the
%> unit interval @f$[0,1]@f$, @f$q_r@f$, are provided by 
%> <code>quadRule2D()</code> and <code>quadRule1D()</code>, respectively.
%>
%> All global variables are @f$\#\mathcal{P} \times 1@f$ <code>cell</code>-arrays, 
%> with @f$\mathcal{P} = \{2p, 2p+1\}@f$ the set of required polynomial orders. 
%> Note, for @f$p=0@f$ only order @f$1@f$ is provided.
%>
%> This function computes the following global variables (dimensions given for
%> each order):
%> - <code>gPhi1D</code>: @f$\hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(q_r)
%>                           \; [R \times N \times 3]@f$
%> - <code>gThetaPhi1D</code>: @f$\hat{\varphi}_i \circ 
%>                                \hat{\mathbf{\gamma}}_{n^-} \circ
%>                                \hat{\vartheta}_{n^-n^+}(q_r)
%>                                \; [R \times N \times 3 \times 3]@f$
%> - <code>gPhi2D</code>: @f$\hat{\varphi}_i(\mathbf{q}_r) \; [R \times N]@f$
%> - <code>gGradPhi2D</code>: @f$\hat{\nabla} \hat{\varphi}_i(\mathbf{q}_r)
%>                               \; [R \times N \times 2]@f$
%> 
%> @param  N          The number of local degrees of freedom. For polynomial
%>                    order @f$p@f$, it is given as @f$N = (p+1)(p+2)/2@f$
%>                    @f$[\text{scalar}]@f$
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
function computeBasesOnQuad(N)
global gPhi2D gGradPhi2D gPhi1D gThetaPhi1D

% Check for valid number of DOFs: N == (p+1)(p+2)/2
assert(~isempty(find(N == ((0:4)+1).*((0:4)+2)/2, 1)), 'Number of degrees of freedom does not match a polynomial order') 

% Determine polynomial degree and quadrature orders
p = (sqrt(8*N+1)-3)/2;
if p > 0, requiredOrders = [2*p, 2*p+1]; else requiredOrders = 1; end

% Initialize global variables
gPhi2D = cell(max(requiredOrders),1);  gGradPhi2D  = cell(max(requiredOrders),1);
gPhi1D = cell(max(requiredOrders),1);  gThetaPhi1D = cell(max(requiredOrders),1);

% Fill global variables
for it = 1 : length(requiredOrders)
  ord = requiredOrders(it);
  [Q1, Q2, ~] = quadRule2D(ord);
  gPhi2D{ord}      = zeros(length(Q1), N);
  for i = 1 : N
    gPhi2D{ord}(:, i) = phi(i, Q1, Q2);
  end % for
  gGradPhi2D{ord}  = zeros(length(Q1), N, 2);
  for m = 1 : 2
    for i = 1 : N
      gGradPhi2D{ord}(:, i, m) = gradPhi(i, m, Q1, Q2);
    end % for
  end % for
  [Q, ~] = quadRule1D(ord);
  gPhi1D{ord} = zeros(length(Q), N, 3);
  for nn = 1 : 3
    [Q1, Q2] = gammaMap(nn, Q);
    for i = 1 : N
      gPhi1D{ord}(:, i, nn) = phi(i, Q1, Q2);
    end
    for np = 1 : 3
      [QP1,QP2] = theta(nn, np, Q1, Q2);
      for i = 1 : N
        gThetaPhi1D{ord}(:, i, nn, np) = phi(i, QP1, QP2);
      end % for
    end % for
  end % for
end % for
end % function
