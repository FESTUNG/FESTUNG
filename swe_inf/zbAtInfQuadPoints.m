
function zbInf = zbAtInfQuadPoints(gInf, zbAlg, p, pInf, beta)
ord    = max(2*p   ,1);
ordInf = max(2*pInf,1);  [Q1, Q2, ~] = quadRuleInf2D(ord, ordInf, beta);
%F: Interpoliert in X2-Richtung linear und ist in X1-Richtung konstant (auf
%Referenzelement); diese Eigenschaft bleibt auf dem physikalischen Element 
%im lokalen Koordinatenystem erhalten => ausreichend auf dem Referenzelem
F = @(X1, X2) (zbAlg(gInf.coordV(gInf.V0Inf(:,2),1),gInf.coordV(gInf.V0Inf(:,2),2)) - ...
  zbAlg(gInf.coordV(gInf.V0Inf(:,1),1),gInf.coordV(gInf.V0Inf(:,1),2)))*X2 ...
  + zbAlg(gInf.coordV(gInf.V0Inf(:,1),1),gInf.coordV(gInf.V0Inf(:,1),2))*ones(size(X2));
zbInf = F(Q1,Q2);
end % function
