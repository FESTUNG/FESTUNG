function ret = assembleVecElemPhiFlux( g, N, sourceEval, basesOnQuad )
%Assert
K = g.numT;
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, ~, W] = quadRule2D(qOrd);
% Assemble matrix
ret = zeros( K*N, 1 );
[~,R] = size(sourceEval);
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
validateattributes(sourceEval, {'numeric'}, {'size', [K R]}, mfilename, 'fluxEval');

%Matlab 2016b
% fluxNu2016b = fluxEval( :, :, 1, :) .* g.nuE0T( :, :, 1 ) + fluxEval( :, :, 2, :) .* g.nuE0T( :, :, 2 );
%Matlab 2016a
% fluxNu = fluxEval( :, :, 1, :) .* repmat( g.nuE0T( :, :, 1 ) , 1, 1, 1, R) ...
%          + fluxEval( :, :, 2, :) .* repmat( g.nuE0T( :, :, 2 ) , 1, 1, 1, R );
     
% ret = W .*  basesOnQuad.phi2D{qOrd}( :, :);
% ret = W' .* basesOnQuad.phi2D{qOrd}( :, :)

ret = reshape( (sourceEval * (W' .* basesOnQuad.phi2D{qOrd}( :, :)))', N*K, 1);
% for n=1:3
%     tmp = zeros( K*N, 1 );
%     for r=1:R
%         tmp(:) = tmp(:) + kron( fluxNu( :, n, 1, r), speye(N) ) * (W(r) .* basesOnQuad.phi1D{qOrd}( r, :, n)' );
%     end
%     fac = g.areaE0T( :, n ) .* markE0Tbdr(:, n);
%     ret(:) = ret(:) + sum( kron( fac, speye(N) ), 2) .* tmp;
% end
end