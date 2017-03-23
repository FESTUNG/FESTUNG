function ret = assembleMatElemPhiDphiFlux( g, N, uEval, Ghat)
%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;
[~, ~, R] = size(Ghat{1});
validateattributes(Ghat, {'cell'}, {'size', [2 1]}, mfilename, 'Ghat');
validateattributes(Ghat{1}, {'numeric'}, {'size', [N N R]}, mfilename, 'Ghat{1}');
validateattributes(Ghat{2}, {'numeric'}, {'size', [N N R]}, mfilename, 'Ghat{2}');
validateattributes(uEval, {'numeric'}, {'size', [K R 2]}, mfilename, 'uEval');

% Assemble matrix
ret = cell(2,1);
ret{1} = sparse(K*N, K*N);
ret{2} = sparse(K*N, K*N);
for r = 1:R
    ret{1} = ret{1} ...
            + kron(spdiags(g.B(:,2,2) .* uEval( :, r, 1), 0,K,K), Ghat{1}(:, :, r) ) ...
            - kron(spdiags(g.B(:,2,1) .* uEval( :, r, 1), 0,K,K), Ghat{2}(:, :, r) );
    ret{2} = ret{2} ...
            - kron(spdiags(g.B(:,1,2) .* uEval( :, r, 2), 0,K,K), Ghat{1}(:, :, r) ) ...
            + kron(spdiags(g.B(:,1,1) .* uEval( :, r, 2), 0,K,K), Ghat{2}(:, :, r) );
end %for
end