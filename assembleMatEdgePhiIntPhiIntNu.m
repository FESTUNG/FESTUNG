function [ globQN ] = assembleGlobQNsupra( g , markE0Tbdr , hatQdiag )

K = g.numTsupra;
N = size(hatQdiag, 1);
globQN = cell(2, 1);
globQN{1} = sparse(K*N, K*N);
globQN{2} = sparse(K*N, K*N);

for n = 1 : 4
    globQN{1} = globQN{1} + kron( spdiags(markE0Tbdr(:,n) .* g.NuLengthE0Tsupra(:,n,1), ...
        0, K, K) , hatQdiag(:,:,n) );
    globQN{2} = globQN{2} + kron( spdiags(markE0Tbdr(:,n) .* g.NuLengthE0Tsupra(:,n,2), ...
        0, K, K) , hatQdiag(:,:,n) );
end  % for

end  % function