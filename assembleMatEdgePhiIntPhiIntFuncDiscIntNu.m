function [ globRD ] = assembleGlobRDsub( g , markE0Tbdr , hatRdiag , K11DG , K12DG , K21DG , K22DG )

[K, N] = size(K11DG);
globRD = cell(2,1);
globRD{1} = sparse(K*N, K*N);
globRD{2} = sparse(K*N, K*N);

for n = 1 : 4
    for i = 1 : N
        globRD{1} = globRD{1} + kron( spdiags( markE0Tbdr(:,n) .* (g.NuLengthE0Tsub(:,n,1) .* K11DG(:,i) ...
            + g.NuLengthE0Tsub(:,n,2) .* K21DG(:,i)) , 0 , K , K ) , hatRdiag(:,:,i,n) );
        globRD{2} = globRD{2} + kron( spdiags( markE0Tbdr(:,n) .* (g.NuLengthE0Tsub(:,n,1) .* K12DG(:,i) ...
            + g.NuLengthE0Tsub(:,n,2) .* K22DG(:,i)) , 0 , K , K ) , hatRdiag(:,:,i,n) );
    end  % for i
end  % for n

end  % function