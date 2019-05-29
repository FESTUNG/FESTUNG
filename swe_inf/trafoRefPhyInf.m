function [X1P, X2P] = trafoRefPhyInf (X1,X2, idEdge, M, H, g, coordAxInf)
ed = g.coordV(g.V0E(idEdge,2),:)-g.coordV(g.V0E(idEdge,1),:);
xHat = g.coordV(g.V0E(idEdge,1),:) + 0.5*ed;
X1PLoc = X1;
X2PLoc = (H(idEdge)+2*X1*M(idEdge)).*(X2-0.5);
x = coordAxInf(idEdge,:);
alpha = (1-2*(x(2)<0))*acos(x(1)/norm(x));%Verdrehung der lok. x-Achse zur globalen
beta = (1-2*(ed(1)>0))*acos(ed(2)/norm(ed));%Verdrehung der lok. y-Achse zur globalen
X1P = xHat(1) + cos(alpha)*X1PLoc - sin(beta)*X2PLoc;
X2P = xHat(2) + sin(alpha)*X1PLoc + cos(beta)*X2PLoc;
end%function