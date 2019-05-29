function detTrafo = trafoDetInf (M, H, g, coordAxInf, markE0TOS)
detTrafo = cell(g.numT,1);
for nT = 1 : g.numT
  for i = 1:3
    idEdge = g.E0T(nT,i);
    if ~markE0TOS(nT,i)
      detTrafo{idEdge,1} = @(x1, x2) zeros(size(x1));
      continue;
    end %if
    ed = g.coordV(g.V0E(idEdge,2),:)-g.coordV(g.V0E(idEdge,1),:);
    m = M(idEdge);
    h = H(idEdge);
    x = coordAxInf(idEdge,:);
    alpha = (1-2*(x(2)<0))*acos(x(1)/norm(x));%Verdrehung der lok. x-Achse zur globalen
    beta = (1-2*(ed(1)>0))*acos(ed(2)/norm(ed));%Verdrehung der lok. y-Achse zur globalen
    detTrafo{idEdge,1} = @(x,y) cos(alpha)*cos(beta)*(h+2*x*m)...
      + sin(beta)*(h+2*x*m).*(sin(alpha)+2*cos(beta)*m*(y-0.5));
  end%for
end%for
end%function