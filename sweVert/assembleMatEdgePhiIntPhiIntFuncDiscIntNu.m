function [ globRD ] = assembleGlobRbdrSupra( g , markE0Tbdr , hatRdiag , uDG )

[K, N] = size(uDG);
globRD = cell(2,1);
globRD{1} = sparse(K*N, K*N);
globRD{2} = sparse(K*N, K*N);

for n = 1 : 4
    for i = 1 : N
        globRD{1} = globRD{1} + kron( spdiags( markE0Tbdr(:,n) .* g.NuLengthE0Tsupra(:,n,1) .* uDG(:,i) ...
            , 0 , K , K ) , hatRdiag(:,:,i,n) );
        globRD{2} = globRD{2} + kron( spdiags( markE0Tbdr(:,n) .* g.NuLengthE0Tsupra(:,n,2) .* uDG(:,i) ...
            , 0 , K , K ) , hatRdiag(:,:,i,n) );
    end  % for i
end  % for n

end  % function