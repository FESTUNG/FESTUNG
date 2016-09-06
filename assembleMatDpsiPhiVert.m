function [ globIntH ] = assembleGlobIntH2( g , quotH, gradPsiPhiC , gradPsiPhiX )

K = g.numTsupra;
N = size(gradPsiPhiC, 1);
R = size(quotH, 2);

globIntH = sparse(K*N, K*size(gradPsiPhiC,2));

for r1 = 1 : R
    helper = repmat(quotH(:,r1), g.NZsupra, 1);
%    helper = repmat(ones(size(quotH,1),1), g.NZsupra, 1);
    for r2 = 1 : R
        globIntH = globIntH - kron( spdiags(g.DAsupra .* helper, 0, K, K) , gradPsiPhiC(:,:,(r1-1)*R+r2) ) ...
            - kron( spdiags(g.ACBDsupra .* helper, 0, K, K) , gradPsiPhiX(:,:,(r1-1)*R+r2) );
    end  % for r2
end  % for r1

globIntH = repmat( speye(g.NX * N), 1 , g.NZsupra ) * globIntH;

end  % function