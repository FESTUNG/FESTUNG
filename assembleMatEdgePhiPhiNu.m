function [ globQ ] = assembleGlobQsub( g , markE0Tint , hatQdiag , hatQoffdiag )

K = g.numTsub;
N = size(hatQdiag, 1);
globQ = cell(2,1);
globQ{1} = sparse(K*N, K*N);
globQ{2} = sparse(K*N, K*N);

for n = 1 : 4
    globQ{1} = globQ{1} + kron( spdiags( 0.5 * markE0Tint(:,n) .* g.NuLengthE0Tsub(:,n,1), ...
        0, K, K ) , hatQdiag(:,:,n) ) + kron( bsxfun(@times, 0.5 * g.markE0TE0Tsub{n}, ...
        g.NuLengthE0Tsub(:,n,1)) , hatQoffdiag(:,:,n) );
    globQ{2} = globQ{2} + kron( spdiags( 0.5 * markE0Tint(:,n) .* g.NuLengthE0Tsub(:,n,2), ...
        0, K, K ) , hatQdiag(:,:,n) ) + kron( bsxfun(@times, 0.5 * g.markE0TE0Tsub{n}, ...
        g.NuLengthE0Tsub(:,n,2)) , hatQoffdiag(:,:,n) );
end  % for

end  % function