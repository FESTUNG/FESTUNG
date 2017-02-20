function ret = assembleVecEdgePhiIntVal( g, N, cEval, markE0Tbdr, basesOnQuad )
%Assert
K = g.numT;
Kedge = g.numE;

validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
R = size( W, 2 );
% Assemble matrix
ret = zeros( K*N, 1 );
ret2016b = zeros( K*N, 1 );

for n=1:3
    fac = g.areaE0T( :, n ) .* markE0Tbdr(:, n);
    % Matlab 2016b
%     ret2016b(:) =ret2016b(:) + reshape( ( fac .* cEval(:,:, n) * (W' .* basesOnQuad.phi1D{qOrd}( :, :, n) ))', K*N, 1 );
    % Matlab 2016a
    ret(:) =ret(:) + reshape( ( repmat(fac, 1, R) .* cEval(:,:, n) * (repmat( W', 1, N).* basesOnQuad.phi1D{qOrd}( :, :, n)) )', K*N, 1 );
end

end
