function [ globIntQ ] = assembleGlobIntQ( g , markE0Tint , hatIntQdiag , hatIntQoffdiag )

K = g.numTsupra;
N1 = size(hatIntQdiag, 1);
N2 = size(hatIntQdiag, 2);
globIntQ = sparse(K*N1, K*N2);

for n = 1 : 2
    globIntQ = globIntQ + kron( spdiags( 0.5 * markE0Tint(:,n+2) .* g.NuLengthE0Tsupra(:,n+2,1), ...
        0, K, K ) , hatIntQdiag(:,:,n) ) + kron( bsxfun(@times, 0.5 * g.markE0TE0Tsupra{n+2}, ...
        g.NuLengthE0Tsupra(:,n+2,1)) , hatIntQoffdiag(:,:,n) );
end  % for

globIntQ = repmat( speye(g.NX * size(hatIntQdiag,1)), 1 , g.NZsupra ) * globIntQ;

end  % function