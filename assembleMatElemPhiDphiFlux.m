function ret = assembleMatElemPhiDphiFlux( g, N, uEval, Gbar)
%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;
[~, ~, R] = size(Gbar{1});
validateattributes(Gbar, {'cell'}, {'size', [2 1]}, mfilename, 'Gbar');
validateattributes(Gbar{1}, {'numeric'}, {'size', [N N R]}, mfilename, 'Gbar{1}');
validateattributes(Gbar{2}, {'numeric'}, {'size', [N N R]}, mfilename, 'Gbar{2}');
validateattributes(uEval, {'numeric'}, {'size', [K R 2]}, mfilename, 'uEval');

% Assemble matrix
ret = cell(2,1);
ret{1} = sparse(K*N, K*N);
ret{2} = sparse(K*N, K*N);
for r = 1:R
    ret{1} = ret{1} ...
            + kron(spdiags(g.B(:,2,2) .* uEval( :, r, 1), 0,K,K), Gbar{1}(:, :, r) ) ...
            - kron(spdiags(g.B(:,2,1) .* uEval( :, r, 1), 0,K,K), Gbar{2}(:, :, r) );
    ret{2} = ret{2} ...
            - kron(spdiags(g.B(:,1,2) .* uEval( :, r, 2), 0,K,K), Gbar{1}(:, :, r) ) ...
            + kron(spdiags(g.B(:,1,1) .* uEval( :, r, 2), 0,K,K), Gbar{2}(:, :, r) );
end %for
end