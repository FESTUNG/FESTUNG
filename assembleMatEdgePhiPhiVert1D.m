function [ globIntW ] = assembleGlobIntW( g , markE0Tint , hatIntW , uDG )

K = g.numTsupra;
N1 = size(hatIntW, 1);
N2 = size(hatIntW, 2);
globIntW = cell(2, 1);
globIntW{1} = sparse(K*N1, K*N2);
globIntW{2} = sparse(K*N1, K*N2);

for i = 1 : size(uDG,2)
    globIntW{1} = globIntW{1} + kron( spdiags(0.5 * markE0Tint(:,1+2) .* g.lengthE0Tsupra(:,1+2) .* uDG(:,i) , ...
        0, K, K) , hatIntW(:,:,i,1) );
    globIntW{2} = globIntW{2} + kron( spdiags(0.5 * markE0Tint(:,2+2) .* g.lengthE0Tsupra(:,2+2) .* uDG(:,i) , ...
        0, K, K) , hatIntW(:,:,i,2) );
end  % for i


globIntW{1} = repmat( speye(g.NX * size(hatIntW,1)), 1 , g.NZsupra ) * globIntW{1};
globIntW{2} = repmat( speye(g.NX * size(hatIntW,1)), 1 , g.NZsupra ) * globIntW{2};

end  % function