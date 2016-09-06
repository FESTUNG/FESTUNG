function [ globM1D ] = assembleGlobM1D( g , hatM1D )

globM1D = kron( g.deltaX * speye(g.NX) , hatM1D );

end  % function