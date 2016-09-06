function [ globIntQN ] = assembleGlobIntQN( g , markE0Tbdr , hatQdiag )

K = g.numTsupra;

globIntQN = kron( spdiags(markE0Tbdr(:,3) .* g.NuLengthE0Tsupra(:,3,1), ...
        0, K, K) , hatQdiag(:,:,1) ) ...
        + kron( spdiags(markE0Tbdr(:,4) .* g.NuLengthE0Tsupra(:,4,1), ...
        0, K, K) , hatQdiag(:,:,2) );
    
globIntQN = repmat( speye(g.NX * size(hatQdiag,1)), 1 , g.NZsupra ) * globIntQN;

end  % function