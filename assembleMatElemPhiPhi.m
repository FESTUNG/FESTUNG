function [ globM ] = assembleGlobMsub( g , hatMc , hatMx )

globM = kron( spdiags(g.deltaX * g.DAsub, 0, g.numTsub, g.numTsub) , hatMc ) ...
    + kron( spdiags(g.deltaX * g.ACBDsub, 0, g.numTsub, g.numTsub) , hatMx );

end  % function