function [ globQTopDown ] = assembleGlobQsupraTopDown( g , markE0Tint , markE0Tfree , hatQdiag , hatQoffdiag )

K = g.numTsupra;
N = size(hatQdiag, 1);
globQTopDown = sparse(K*N, K*N);

for n = 1 : 2
    globQTopDown = globQTopDown + kron( spdiags( 0.5 * markE0Tint(:,n) .* g.NuLengthE0Tsupra(:,n,2), ...
        0, K, K ) , hatQdiag(:,:,n) ) + kron( bsxfun(@times, 0.5 * g.markE0TE0Tsupra{n}, ...
        g.NuLengthE0Tsupra(:,n,2)) , hatQoffdiag(:,:,n) );
end  % for

globQTopDown = globQTopDown + kron( spdiags( markE0Tfree(:,2) .* g.NuLengthE0Tsupra(:,2,2), ...
     0, K, K) , hatQdiag(:,:,2) );
 
end  % function