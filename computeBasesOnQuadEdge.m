% TODO
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
