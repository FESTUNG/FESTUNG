function [ globCbdr ] = assembleGlobCbdrSupra( g , markE0Tbdr , hatCdiag )

K = g.numTsupra;
N1 = size(hatCdiag, 1);
N2 = size(hatCdiag, 2);
globCbdr = cell(2, 1);
globCbdr{1} = sparse(K*N1, K*N2);
globCbdr{2} = sparse(K*N1, K*N2);

for n = 1 : 4
    globCbdr{1} = globCbdr{1} + kron( spdiags(markE0Tbdr(:,n) .* g.NuLengthE0Tsupra(:,n,1), ...
        0, K, K) , hatCdiag(:,:,n) );
    globCbdr{2} = globCbdr{2} + kron( spdiags(markE0Tbdr(:,n) .* g.NuLengthE0Tsupra(:,n,2), ...
        0, K, K) , hatCdiag(:,:,n) );
end  % for

end  % function