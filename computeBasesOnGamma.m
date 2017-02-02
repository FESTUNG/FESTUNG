%
function basesOnGamma = computeBasesOnGamma(N, basesOnGamma, requiredOrders)
% Check for valid number of DOFs: N == (p+1)(p+2)/2
assert(~isempty(find(N == (0:3)+1, 1)), ...
  'Number of degrees of freedom does not match a polynomial order')
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma')

% Determine polynomial degree and quadrature orders
% if nargin < 3
%   p = (sqrt(8*N+1)-3)/2;
%   if p > 0
%     requiredOrders = [2*p, 2*p+1]; 
%   else
%     requiredOrders = 1; 
%   end % if
% end % if

if nargin < 3
  p = N - 1;
  if p > 0
    requiredOrders = [2*p, 2*p+1]; 
  else
    requiredOrders = 1; 
  end % if
end % if

% Initialize global variables
basesOnGamma.phi1D = cell(max(requiredOrders),1);  
basesOnGamma.thetaPhi1D = cell(max(requiredOrders),1);
% Fill global variables
for it = 1 : length(requiredOrders)
  ord = requiredOrders(it);

  [Q, ~] = quadRule1D(ord);
  basesOnGamma.phi1D{ord} = zeros(length(Q), N, 3);
  for nn = 1 : 3
    [Q1, Q2] = gammaMap(nn, Q);
    for i = 1 : N
      basesOnGamma.phi1D{ord}(:, i, nn) = phi1D(i, Q);
    end
%     for np = 1 : 3
%       [QP1,QP2] = theta(nn, np, Q1, Q2);
%       for i = 1 : N
%         basesOnGamma.thetaPhi1D{ord}(:, i, nn, np) = phi(i, QP1, QP2);
%       end % for
%     end % for
  end % for
end % for


% % Fill global variables
% for it = 1 : length(requiredOrders)
%   ord = requiredOrders(it);
%   [Q1, Q2, ~] = quadRule2D(ord);
%   basesOnGamma.phi2D{ord}      = zeros(length(Q1), N);
%   for i = 1 : N
%     basesOnGamma.phi2D{ord}(:, i) = phi(i, Q1, Q2);
%   end % for
%   basesOnGamma.gradPhi2D{ord}  = zeros(length(Q1), N, 2);
%   for m = 1 : 2
%     for i = 1 : N
%       basesOnGamma.gradPhi2D{ord}(:, i, m) = gradPhi(i, m, Q1, Q2);
%     end % for
%   end % for
%   [Q, ~] = quadRule1D(ord);
%   basesOnGamma.phi1D{ord} = zeros(length(Q), N, 3);
%   for nn = 1 : 3
%     [Q1, Q2] = gammaMap(nn, Q);
%     for i = 1 : N
%       basesOnGamma.phi1D{ord}(:, i, nn) = phi(i, Q1, Q2);
%     end
%     for np = 1 : 3
%       [QP1,QP2] = theta(nn, np, Q1, Q2);
%       for i = 1 : N
%         basesOnGamma.thetaPhi1D{ord}(:, i, nn, np) = phi(i, QP1, QP2);
%       end % for
%     end % for
%   end % for
% end % for
end % function
