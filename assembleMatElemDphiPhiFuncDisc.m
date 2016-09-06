function [ globG ] = assembleGlobGsupra( g , hatGc , hatGx , hatGy , uDG )

[K, N] = size(uDG);
globG = cell(2,1);
globG{1} = sparse(K*N, K*N);
globG{2} = sparse(K*N, K*N);

for i = 1 : N
    globG{1} = globG{1} - kron(spdiags(uDG(:,i) .* g.ACBDsupra, 0, K, K), hatGx(:,:,i,1)) ...
                        - kron(spdiags(uDG(:,i) .* g.DAsupra  , 0, K, K), hatGc(:,:,i,1)) ...
                        + kron(spdiags(uDG(:,i) .* g.ACBDsupra, 0, K, K), hatGy(:,:,i,2)) ...
                        + kron(spdiags(uDG(:,i) .* g.BAsupra  , 0, K, K), hatGc(:,:,i,2));
    globG{2} = globG{2} - kron(spdiags(uDG(:,i) .* g.deltaX   , 0, K, K), hatGc(:,:,i,2));
end  % for

end  % function