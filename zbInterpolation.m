function zbAlg = zbInterpolation(g, depth, X1, X2)
%% gives linear continuous interpolation of zbAlg
lenX = size(X1,2);
depth = depth(:);
FinvX1 = 0.5 ./ (g.areaT * ones(1,lenX)) .* ( bsxfun(@times, g.B(:,2,2), X1) - bsxfun(@times, g.B(:,1,2), X2) - (g.B(:,2,2).*g.coordV0T(:,1,1) - g.B(:,1,2).*g.coordV0T(:,1,2) ) * ones(1,lenX) );
FinvX2 = 0.5 ./ (g.areaT * ones(1,lenX)) .* (-bsxfun(@times, g.B(:,2,1), X1) + bsxfun(@times, g.B(:,1,1), X2) + (g.B(:,2,1).*g.coordV0T(:,1,1) - g.B(:,1,1).*g.coordV0T(:,1,2) ) * ones(1,lenX) );
zbAlg  = bsxfun(@times, depth(g.V0T(:,1)), 1 - FinvX1 - FinvX2) ...
       + bsxfun(@times, depth(g.V0T(:,2)),     FinvX1         ) ...
       + bsxfun(@times, depth(g.V0T(:,3)),             FinvX2 );
end % function