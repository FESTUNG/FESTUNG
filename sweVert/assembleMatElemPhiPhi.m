function [ globM ] = assembleGlobMsupra( g , hatMc , hatMx )

globM = kron( spdiags(g.deltaX * g.DAsupra, 0, g.numTsupra, g.numTsupra) , hatMc ) ...
    + kron( spdiags(g.deltaX * g.ACBDsupra, 0, g.numTsupra, g.numTsupra) , hatMx );

end  % function