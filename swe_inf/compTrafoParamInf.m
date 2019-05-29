function [M, H, coordAxInf, Alpha, Beta, xHat] = compTrafoParamInf (gInf)
%Problem: alle Referenzpunkte brauchen in Normalenrichtung den gleichen
%Abstand damit die Transformationsformel hergeleitet werden kann!
M = zeros(gInf.numInfElem,1);
H = zeros(gInf.numInfElem,1);
Alpha = zeros(gInf.numInfElem,1);
Beta = zeros(gInf.numInfElem,1);
coordAxInf = zeros(gInf.numInfElem,2);
xHat = zeros(gInf.numInfElem,2);
for nT = 1: gInf.numInfElem
  a = gInf.infDir(gInf.V0Inf(nT,2),:);
  b = gInf.infDir(gInf.V0Inf(nT,1),:);
  nu(1) = -gInf.nuE0Inf(nT,2,1);
  nu(2) = -gInf.nuE0Inf(nT,2,2);
  ed = gInf.coordV(gInf.V0Inf(nT,2),:)-gInf.coordV(gInf.V0Inf(nT,1),:);
  A=[nu',ed'];

  coeff = A\a';
  a = a/coeff(1);
  coeff = A\b';
  b=b/coeff(1);
  x = 0.5*(a+b);
  xHat(nT,:) = gInf.coordV(gInf.V0Inf(nT,1),:) + 0.5*ed;
  if a == b
    M(nT)= 0;
  else      
    xTilde = xHat(nT, :)+x;
    A = [a',-x'];
    coeff = A\(xTilde-gInf.coordV(gInf.V0Inf(nT,2),:))';
    xBar = gInf.coordV(gInf.V0Inf(nT,2),:) + coeff(1)*a;
    M(nT) = norm(gInf.coordV(gInf.V0Inf(nT,2),:) + a -xTilde)/norm(xTilde-xBar);
  end%if
  H(nT) = norm(ed);
  coordAxInf(nT,1) = x(1);
  coordAxInf(nT,2) = x(2);
  coordAxInf(nT,:) = coordAxInf(nT,:)/norm(coordAxInf(nT,:));
  Alpha(nT) = (1-2*(x(2)<0))*acos(x(1)/norm(x));%Verdrehung der lok. x-Achse zur globalen
  Beta(nT) = (1-2*(ed(1)>0))*acos(ed(2)/norm(ed));%Verdrehung der lok. y-Achse zur globalen
  
end % for