function [ globR ] = assembleGlobRsub( g , markE0Tint , hatRdiag , hatRoffdiag , ...
    K11DG , K12DG , K21DG , K22DG )

[K, N] = size(K11DG);
globR = cell(2, 1);
globR{1} = sparse(K*N, K*N);
globR{2} = sparse(K*N, K*N);

for n = 1 : 4
    for i = 1 : N
        helper1 = 0.5 * g.NuLengthE0Tsub(:,n,1) * K11DG(:,i)' + 0.5 * g.NuLengthE0Tsub(:,n,2) * K21DG(:,i)';
        helper2 = 0.5 * g.NuLengthE0Tsub(:,n,1) * K12DG(:,i)' + 0.5 * g.NuLengthE0Tsub(:,n,2) * K22DG(:,i)';
        globR{1} = globR{1} + kron( g.markE0TE0Tsub{n} .* helper1 , hatRoffdiag(:,:,i,n) ) ...
            + kron( spdiags( markE0Tint(:,n) .* diag(helper1) , 0, K, K) , hatRdiag(:,:,i,n) );
        globR{2} = globR{2} + kron( g.markE0TE0Tsub{n} .* helper2 , hatRoffdiag(:,:,i,n) ) ...
            + kron( spdiags( markE0Tint(:,n) .* diag(helper2) , 0, K, K) , hatRdiag(:,:,i,n) );
    end  % for i
end  % for n

end  % function