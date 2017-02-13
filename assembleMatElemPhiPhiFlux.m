function ret = assembleMatElemPhiPhiFlux( g, N, uEval, Gbar )
%Assert that the length of uEval and hatSbarOnQuad is ok
K = g.numT;

% Assemble matrix
ret = cell(2,1); 
ret{1} = sparse(K*N, K*N); 
ret{2} = sparse(K*N, K*N);

% warning('flipping in assembleMatElemPhiPhiFlux');

% ret( i, j, ip, m )

for iT = 1:K
    iTs = (iT-1)*N + 1;
    iTe = (iT)*N;
    
    for iDim = 1:2
        tmp = zeros(N, N);
        for i = 1:N
            for j=1:N
%                 Gbar( :, i, j, iDim)
%                 uEval( iT, :, iDim)
                tmp(i,j) = uEval( iT, :, iDim) * Gbar( :, i, j, iDim) ;
%                 tmp(i,j) = fliplr(uEval( iT, :, iDim)) * Gbar( :, i, j, iDim) ;
            end
        end
        ;
        ret{iDim}( iTs:iTe,  iTs:iTe) = ret{iDim}( iTs:iTe,  iTs:iTe ) + 2 .* g.areaT( iT ) .* tmp;
    end
end

end