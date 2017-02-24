%
function basesOnGamma = computeBasesOnGamma(N, basesOnGamma, requiredOrders)
% Check for valid number of DOFs: N == (p+1)(p+2)/2
assert(~isempty(find(N == (0:4)+1, 1)), ...
    'Number of degrees of freedom does not match a polynomial order')
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma')

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
    basesOnGamma.phi1D{ord} = zeros(length(Q), N);
    for i = 1 : N
        basesOnGamma.phi1D{ord}(:, i) = phi1D(i, Q);
    end
    
    basesOnGamma.thetaPhi1D{ord} = zeros(length(Q), N, 2);
    basesOnGamma.thetaPhi1D{ord}(:, :, 1) = basesOnGamma.phi1D{ord}(:, :);
    basesOnGamma.thetaPhi1D{ord}(:, :, 2) = flipud(basesOnGamma.phi1D{ord}(:, :));
end % for
end % function
