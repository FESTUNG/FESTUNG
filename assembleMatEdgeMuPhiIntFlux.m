function ret = assembleMatEdgeMuPhiIntFlux( g, markE0Tbdr, uEdge, Sbar )
%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;
KEdge = g.numE;
[N, Nmu, ~, R] = size(Sbar{1});
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(Sbar, {'cell'}, {'size', [2 1]}, mfilename, 'Sbar');
validateattributes(Sbar{1}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'Sbar{1}');
validateattributes(Sbar{2}, {'numeric'}, {'size', [N Nmu 3 R]}, mfilename, 'Sbar{2}');
validateattributes(uEdge, {'cell'}, {'size', [2 1]}, mfilename, 'uEdge');
validateattributes(uEdge{1}, {'numeric'}, {'size', [R K 3]}, mfilename, 'uEdge{1}');
validateattributes(uEdge{2}, {'numeric'}, {'size', [R K 3]}, mfilename, 'uEdge{2}');

% Assemble matrix
ret = cell(2,1);
ret{1} = sparse( K*N, KEdge*Nmu);
ret{2} = sparse( K*N, KEdge*Nmu);
for iDim = 1:2
    for iE = 1:3
        for l=1:2
        Rkn = markE0Tbdr(:,iE) .* g.areaE0T(:, iE) .*  g.markSideE0T(:, iE, l) .* g.nuE0T( :, iE , iDim);
            for r = 1:R
                ret{iDim} = ret{iDim} + ...
                    kron( sparse( 1:K, g.E0T(:, iE), ones(K, 1) .* uEdge{iDim}( r, :, iE)' .* Rkn , K, KEdge ), ...
                    Sbar{l}( :, :, iE, r) );
            end
        end
    end
end
end
