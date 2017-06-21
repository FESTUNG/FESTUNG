function ret = assembleMatEdgeMuPhiIntFlux( g, markE0Tbdr, uEdge, Shat )
%Assert that the length of uEval and hatShatOnQuad is ok
K = g.numT;
KEdge = g.numE;
[N, Nmu, ~, R] = size(Shat{1});
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(Shat, {'cell'}, {'size', [2 1]}, mfilename, 'Shat');
validateattributes(Shat{1}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'Shat{1}');
validateattributes(Shat{2}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'Shat{2}');
validateattributes(uEdge, {'cell'}, {'size', [2 1]}, mfilename, 'uEdge');
validateattributes(uEdge{1}, {'numeric'}, {'size', [R K 3]}, mfilename, 'uEdge{1}');
validateattributes(uEdge{2}, {'numeric'}, {'size', [R K 3]}, mfilename, 'uEdge{2}');

% Assemble matrix
ret = sparse( K*N, KEdge*Nmu);
for iDim = 1:2
    for iE = 1:3
        for l=1:2
            Rkn = markE0Tbdr(:,iE) .* g.areaE0T(:, iE) .*  g.markSideE0T(:, iE, l) .* g.nuE0T( :, iE , iDim);
            for r = 1:R
                ret = ret + ...
                    kron( sparse( 1:K, g.E0T(:, iE), ones(K, 1) .* uEdge{iDim}( r, :, iE)' .* Rkn , K, KEdge ), ...
                    Shat{l}( :, :, iE, r) );
            end
        end
    end
end
end
