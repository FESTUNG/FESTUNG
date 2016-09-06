function [ globSD ] = assembleGlobSDsub( g , markE0Tbdr , hatSdiag , eta )

K = g.numTsub;
N = size(hatSdiag, 1);
globSD  = sparse(K*N, K*N);

for n = 1 : 4
    globSD = globSD + eta * kron( spdiags(markE0Tbdr(:,n), 0, K, K) , hatSdiag(:,:,n) );
end  % for

end  % function