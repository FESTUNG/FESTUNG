
function zbInfE0I = zbAtInfQuadPointsE0I(gInf, zbAlg, p, beta)
ord = max(2*p+1,1);  [Q1,~] = quadRuleInf1D(ord, beta, 2);
zbInfE0I = zeros(gInf.numInfElem,length(Q1),3);
F = @(X1, i) zbAlg(gInf.coordV(gInf.V0Inf(:,i),1),gInf.coordV(gInf.V0Inf(:,i),2))*ones(size(X1));

zbInfE0I(:,:,1) = F(Q1,2);
zbInfE0I(:,:,3) = F(Q1,1);

end % function
