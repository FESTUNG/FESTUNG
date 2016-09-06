function [ globR ] = assembleGlobRsupra( g , markE0Tint , hatRdiag , hatRoffdiag , uDG )

[K, N] = size(uDG);
globR = cell(2, 1);
globR{1} = sparse(K*N, K*N);
globR{2} = sparse(K*N, K*N);

for n = 1 : 4
    for i = 1 : N
        helper1 = 0.5 * g.NuLengthE0Tsupra(:,n,1) * uDG(:,i)';
        helper2 = 0.5 * g.NuLengthE0Tsupra(:,n,2) * uDG(:,i)';
        globR{1} = globR{1} + kron( g.markE0TE0Tsupra{n} .* helper1 , hatRoffdiag(:,:,i,n) ) ...
            + kron( spdiags( markE0Tint(:,n) .* diag(helper1) , 0, K, K) , hatRdiag(:,:,i,n) );
        globR{2} = globR{2} + kron( g.markE0TE0Tsupra{n} .* helper2 , hatRoffdiag(:,:,i,n) ) ...
            + kron( spdiags( markE0Tint(:,n) .* diag(helper2) , 0, K, K) , hatRdiag(:,:,i,n) );
    end  % for i
end  % for n

end  % function