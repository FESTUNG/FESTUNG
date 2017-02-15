%
function basesOnGamma = computeBasesOnGamma(N, basesOnGamma, requiredOrders)
% Check for valid number of DOFs: N == (p+1)(p+2)/2
assert(~isempty(find(N == (0:3)+1, 1)), ...
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

% warning('flipping bases on Gamma');

% Initialize global variables
basesOnGamma.phi1D = cell(max(requiredOrders),1);
% basesOnGamma.thetaPhi1D = cell(max(requiredOrders),1);
% Fill global variables
for it = 1 : length(requiredOrders)
    ord = requiredOrders(it);
    
    [Q, ~] = quadRule1D(ord);
    basesOnGamma.phi1D{ord} = zeros(length(Q), N);
    for i = 1 : N
        basesOnGamma.phi1D{ord}(:, i) = phi1D(i, Q);
    end
    
    
%     basesOnGamma.thetaPhi1D{ord}(:, i, 1, 1) = basesOnGamma.phi1D{ord}(:, :);
%     basesOnGamma.thetaPhi1D{ord}(:, i, 1, 2) = basesOnGamma.phi1D{ord}(:, :);
%     basesOnGamma.thetaPhi1D{ord}(:, i, 1, 3) = basesOnGamma.phi1D{ord}(:, :);
%     basesOnGamma.thetaPhi1D{ord}(:, i, 2, 1) = fliplr(basesOnGamma.phi1D{ord}(:, :));
%     basesOnGamma.thetaPhi1D{ord}(:, i, 2, 2) = basesOnGamma.phi1D{ord}(:, :);
%     basesOnGamma.thetaPhi1D{ord}(:, i, 2, 3) = basesOnGamma.phi1D{ord}(:, :);
%     basesOnGamma.thetaPhi1D{ord}(:, i, 3, 1) = fliplr(basesOnGamma.phi1D{ord}(:, :));
%     basesOnGamma.thetaPhi1D{ord}(:, i, 3, 2) = fliplr(basesOnGamma.phi1D{ord}(:, :));
%     basesOnGamma.thetaPhi1D{ord}(:, i, 3, 3) = basesOnGamma.phi1D{ord}(:, :);
%     basesOnGamma.phi1D{ord} = fliplr( basesOnGamma.phi1D{ord} );
%     for nn = 1 : 3
%         [Q1, Q2] = gammaMap(nn, Q);
%         for np = 1 : 3
%             [QP1,QP2] = theta(nn, np, Q1, zeros( 1, size(Q1,2)) );
%             for i = 1 : N
%                 %             basesOnGamma.thetaPhi1D{ord}(:, i, nn, np) = phi(i, QP1, QP2);
%                 basesOnGamma.thetaPhi1D{ord}(:, i, nn, np) = phi1D(i, QP1);
%             end % for
%         end % for
%     end
end % for
end % function
