function [ globR1D ] = assembleGlobR1D( g , hatR1D )

K = g.NX;
N = size(hatR1D, 1);
globR1D = cell(2, 1);
globR1D{1} = sparse(K*N, K*N);
globR1D{2} = sparse(K*N, K*N);

globR1D{1} = kron( spdiags( 0.5 * g.markE0Tint1D(:,1), 0, K, K ) , hatR1D(:,:,1) );
globR1D{2} = kron( spdiags( 0.5 * g.markE0Tint1D(:,2), 0, K, K ) , hatR1D(:,:,2) );

end  % function