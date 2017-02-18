function ret = assembleVecEdgePhiIntFlux( g, N, fluxEval, markE0Tbdr, basesOnQuad )
%Assert
K = g.numT;

validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

% Assemble matrix
ret = zeros( K*N, 1 );
[~,~,~,R] = size(fluxEval);

fluxNu = fluxEval( :, :, 1, :) .* g.nuE0T( :, :, 1 ) + fluxEval( :, :, 2, :) .* g.nuE0T( :, :, 2 );

for n=1:3
    tmp = zeros( K*N, 1 );
    for r=1:R
        tmp(:) = tmp(:) + kron( fluxNu( :, n, 1, r), speye(N) ) * (W(r) .* basesOnQuad.phi1D{qOrd}( r, :, n)' );
    end
    fac = g.areaE0T( :, n ) .* markE0Tbdr(:, n);
%     ret(:) = ret(:) + kron( fac .* tmp, speye(N) ) * (W * basesOnQuad.phi1D{qOrd}( :, :, n))';
    ret(:) = ret(:) + sum( kron( fac, speye(N) ), 2) .* tmp;
end

% for n=1:3
%     tmp = zeros( K, 1 );
%     for r=1:R
%         tmp(:) = tmp(:) + fluxEval( :, n, 1, r).* g.nuE0T( :, n, 1 );
%         tmp(:) = tmp(:) + fluxEval( :, n, 2, r).* g.nuE0T( :, n, 2 );
%     end
%     fac = g.areaE0T( :, n ) .* markE0Tbdr(:, n);
%     ret(:) = ret(:) + kron( fac .* tmp, speye(N) ) * (W * basesOnQuad.phi1D{qOrd}( :, :, n))';
% end
end