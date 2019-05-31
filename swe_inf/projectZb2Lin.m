%%WICHTIG: Bei zbLin(:,2) wird NICHT durch H geteilt (vgl. Herleitung von
%%zb(F(xHat))=z_1+(z2-z1)*xHat_2)

function zbLin = projectZb2Lin(gInf, zbAlg)
zbLin = zeros(gInf.numInfElem,2);
zbLin(:,1) = zbAlg(gInf.coordV(gInf.V0Inf(:,1),1),gInf.coordV(gInf.V0Inf(:,1),2));
zbLin(:,2) = (zbAlg(gInf.coordV(gInf.V0Inf(:,2),1),gInf.coordV(gInf.V0Inf(:,2),2)) - ...
  zbAlg(gInf.coordV(gInf.V0Inf(:,1),1),gInf.coordV(gInf.V0Inf(:,1),2)));%./gInf.H;
end % function
