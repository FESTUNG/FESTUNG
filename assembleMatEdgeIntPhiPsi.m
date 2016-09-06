function [ globC ] = assembleGlobCsupra( g , markE0Tint , hatCdiag , hatCoffdiag )

K = g.numTsupra;
N1 = size(hatCdiag, 1);
N2 = size(hatCdiag, 2);
globC = cell(2,1);
globC{1} = sparse(K*N1, K*N2);
globC{2} = sparse(K*N1, K*N2);

for n = 1 : 4
    globC{1} = globC{1} + kron( spdiags( 0.5 * markE0Tint(:,n) .* g.NuLengthE0Tsupra(:,n,1), ...
        0, K, K ) , hatCdiag(:,:,n) ) + kron( bsxfun(@times, 0.5 * g.markE0TE0Tsupra{n}, ...
        g.NuLengthE0Tsupra(:,n,1)) , hatCoffdiag(:,:,n) );
    globC{2} = globC{2} + kron( spdiags( 0.5 * markE0Tint(:,n) .* g.NuLengthE0Tsupra(:,n,2), ...
        0, K, K ) , hatCdiag(:,:,n) ) + kron( bsxfun(@times, 0.5 * g.markE0TE0Tsupra{n}, ...
        g.NuLengthE0Tsupra(:,n,2)) , hatCoffdiag(:,:,n) );
end  % for

end  % function