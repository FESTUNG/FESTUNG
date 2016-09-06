function [ globS1D ] = assembleGlobS1D( g , lambdaMax1D , hatQdiag1D , hatQoffdiag1D )

K = g.NX;
N = size(hatQdiag1D, 1);
globS1D = cell(2,1);
globS1D{1} = sparse(K*N, K*N);
globS1D{2} = sparse(K*N, K*N);

globS1D{1} = kron( spdiags( 0.5 * g.markE0Tint1D(:,1) .* [lambdaMax1D; 0], 0, K, K ) , hatQdiag1D(:,:,1) ) ...
    - kron( 0.5 * spdiags([lambdaMax1D; 0], 0, K, K) * g.markE0TE0T1D{1}, hatQoffdiag1D(:,:,1) );
globS1D{2} = - kron( spdiags( 0.5 * g.markE0Tint1D(:,2) .* [0; lambdaMax1D], 0, K, K ) , hatQdiag1D(:,:,2) ) ...
    + kron( 0.5 * spdiags([0; lambdaMax1D], 0, K, K) * g.markE0TE0T1D{2} , hatQoffdiag1D(:,:,2) );

globS1D = globS1D{1} + globS1D{2};

end  % function