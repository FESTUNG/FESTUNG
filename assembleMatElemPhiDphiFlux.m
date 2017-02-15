function ret = assembleMatElemPhiDphiFlux( g, N, uEval, Gbar )
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

    tmp = zeros(N, N, 2);
    for iDim = 1:2
        for i = 1:N
            for j=1:N
                %                 Gbar( :, i, j, iDim)
                %                 uEval( iT, :, iDim)
                %                 if iT==1
                %                     tmp(i,j, iDim) = uEval( iT, :, iDim) * Gbar( :, i, j, iDim) ;
                %                 else
                %                     tmp(i,j, iDim) = fliplr(uEval( iT, :, iDim)) * Gbar( :, i, j, iDim) ;
                %                 end
                tmp(i,j, iDim) = uEval( iT, :, iDim) * Gbar( :, i, j, iDim) ;
                
                %                 tmp(i,j) = fliplr(uEval( iT, :, iDim)) * Gbar( :, i, j, iDim) ;
            end
        end
        
        
        %         ret{1} = + kron(spdiags(g.B(:,2,2), 0,K,K), refElemDphiPhi(:,:,1)) ...
        %             - kron(spdiags(g.B(:,2,1), 0,K,K), refElemDphiPhi(:,:,2));
        %         ret{2} = - kron(spdiags(g.B(:,1,2), 0,K,K), refElemDphiPhi(:,:,1)) ...
        %             + kron(spdiags(g.B(:,1,1), 0,K,K), refElemDphiPhi(:,:,2));
        
        %         ret{iDim}( iTs:iTe,  iTs:iTe) = ret{iDim}( iTs:iTe,  iTs:iTe ) + 2 .* g.areaT( iT ) .* tmp;
    end
    %     if (iDim == 1)
    %         ret{iDim}( iTs:iTe,  iTs:iTe) = ret{iDim}( iTs:iTe,  iTs:iTe ) + 2 .* tmp;
    ret{1}( iTs:iTe,  iTs:iTe) = ret{1}( iTs:iTe,  iTs:iTe ) + g.B(iT,2,2) .* tmp(:,:, 1) -  g.B(iT,2,1) .* tmp(:,:, 2) ;
    %     else
    ret{2}( iTs:iTe,  iTs:iTe) = ret{2}( iTs:iTe,  iTs:iTe ) - g.B(iT,1,2) .* tmp(:,:, 1) +  g.B(iT,1,1) .* tmp(:,:, 2) ;
    %     end
end

% ret{1}
% ret{2}
end