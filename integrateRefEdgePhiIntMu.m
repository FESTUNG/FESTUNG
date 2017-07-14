% TODO
function ret = integrateRefEdgePhiIntMu(N, basesOnQuad, qOrd)
validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N(1)+1)-3)/2;  qOrd = 2*p+1;  end
[~, W] = quadRule1D(qOrd);

ret = zeros(N(1), N(2), 3, 2); % [N(1) x N(2) x 3 x 2]
for n = 1 : 3
    for i = 1 : N(1)
        for j = 1 : N(2)
            ret(i, j, n, 1) =  W * (basesOnQuad.phi1D{qOrd}(:, i, n) .* basesOnQuad.mu{qOrd}(:, j) );
            ret(i, j, n, 2) =  W * (basesOnQuad.phi1D{qOrd}(:, i, n) .* basesOnQuad.thetaMu{qOrd}(:, j, 2) );
        end % for
    end % for
end % for
end % function
