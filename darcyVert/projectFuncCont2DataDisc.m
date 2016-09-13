function [ fDG ] = projectFuncCont2DataDisc( g , fAlg , p , ord , hatM , basesOnQuad )

[Q1, Q2, W] = gaussQuadRule2D(ord);
F1 = @(X,Y) g.deltaX * ones(g.numTsub,1) * X + g.coordV0Tsub(:,1,1) * ones(size(X));
F2 = @(X,Y) ( g.coordV0Tsub(:,2,2) - g.coordV0Tsub(:,1,2) ) * X ...
    + ( g.coordV0Tsub(:,4,2) - g.coordV0Tsub(:,1,2) ) * Y ...
    + ( g.coordV0Tsub(:,1,2) + g.coordV0Tsub(:,3,2) - g.coordV0Tsub(:,2,2) - g.coordV0Tsub(:,4,2) ) ...
    * (X .* Y) + g.coordV0Tsub(:,1,2) * ones(size(X));

rhs = fAlg( F1(Q1.', Q2.') , F2(Q1.', Q2.') ) * (repmat(W', 1, (p+1)^2) .* basesOnQuad.gPhi2D) ;
fDG = rhs / hatM;

end  % function

