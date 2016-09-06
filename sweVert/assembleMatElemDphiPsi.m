function [ globE ] = assembleGlobEsupra( g , HatEc , HatEx , HatEy )

K = g.numTsupra;
N1 = size(HatEc, 1);
N2 = size(HatEc, 2);
globE = cell(2, 1);
globE{1} = sparse(K*N1, K*N2);
globE{2} = sparse(K*N1, K*N2);

globE{1} =   kron( spdiags(g.BAsupra  , 0, K, K) , HatEc(:,:,2) ) ...
           + kron( spdiags(g.ACBDsupra, 0, K, K) , HatEy(:,:,2) ) ...
           - kron( spdiags(g.DAsupra  , 0, K, K) , HatEc(:,:,1) ) ...
           - kron( spdiags(g.ACBDsupra, 0, K, K) , HatEx(:,:,1) );
globE{2} = - kron( spdiags(g.deltaX*ones(K,1), 0, K, K) , HatEc(:,:,2) );

% globE{1} = sparse(K*N1, K*N2);

end  % function