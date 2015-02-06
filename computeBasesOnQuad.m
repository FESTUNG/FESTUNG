% This file is part of FESTUNG 
% Copyright (C) 2014 Florian Frank, Balthasar Reuter, Vadym Aizinger
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function computeBasesOnQuad(N)
global gPhi2D gGradPhi2D gPhi1D gThetaPhi1D
p = (sqrt(8*N+1)-3)/2;
if p > 0, requiredOrders = [2*p, 2*p+1]; else requiredOrders = 1; end
gPhi2D = cell(max(requiredOrders),1);  gGradPhi2D  = cell(max(requiredOrders),1);
gPhi1D = cell(max(requiredOrders),1);  gThetaPhi1D = cell(max(requiredOrders),1);
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
