function [ globH ] = assembleGlobHsupra( g , hatHc , hatHx , hatHy )

K = g.numTsupra; N = size(hatHc, 1);
globH = cell(2, 1);
globH{1} = sparse(K*N, K*N);
globH{2} = sparse(K*N, K*N);

globH{1} =   kron( spdiags(g.BAsupra  , 0, K, K) , hatHc(:,:,2) ) ...
           + kron( spdiags(g.ACBDsupra, 0, K, K) , hatHy(:,:,2) ) ...
           - kron( spdiags(g.DAsupra  , 0, K, K) , hatHc(:,:,1) ) ...
           - kron( spdiags(g.ACBDsupra, 0, K, K) , hatHx(:,:,1) );
globH{2} = - kron( spdiags(g.deltaX*ones(K,1), 0, K, K) , hatHc(:,:,2) );

end  % function