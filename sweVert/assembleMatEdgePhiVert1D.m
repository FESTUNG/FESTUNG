function [ globIntV ] = assembleGlobIntV( g , markE0Tint , hatIntV )

K = g.numTsupra;
N1 = size(hatIntV, 1);
N2 = size(hatIntV, 2);
globIntV = cell(2,1);
globIntV{1} = sparse(K*N1, K*N2);
globIntV{2} = sparse(K*N1, K*N2);


globIntV{1} = kron( spdiags( markE0Tint(:,1+2) .* g.NuLengthE0Tsupra(:,1+2), ...
    0, K, K ) , hatIntV(:,:,1) );
globIntV{2} = kron( spdiags( markE0Tint(:,2+2) .* g.NuLengthE0Tsupra(:,2+2), ...
    0, K, K ) , hatIntV(:,:,2) );

globIntV{1} = repmat( speye(g.NX * size(hatIntV,1)), 1 , g.NZsupra ) * globIntV{1};
globIntV{2} = repmat( speye(g.NX * size(hatIntV,1)), 1 , g.NZsupra ) * globIntV{2};

end  % function