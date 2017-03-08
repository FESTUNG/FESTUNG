function ret = assembleVecElemPhiSource( g, N, sourceEval, basesOnQuad )
K = g.numT;
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, ~, W] = quadRule2D(qOrd);
[~,R] = size(sourceEval); %needed for assertion
%Assert
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
validateattributes(sourceEval, {'numeric'}, {'size', [K R]}, mfilename, 'sourceEval');

% Assemble matrix
ret = reshape( (2 * g.areaT .* sourceEval * (W' .* basesOnQuad.phi2D{qOrd}( :, :)))', N*K, 1);
end %function