function infDir = infBdrs(markE0TOS, markE0TL,g)
infDir = zeros(g.numV,2);


for i = 1:3
  for nT = 1:g.numT
    if markE0TOS(nT,i)
      p1 = g.V0E(g.E0T(nT,i),1);
      p2 = g.V0E(g.E0T(nT,i),2);
      norm11 = g.nuE0T(nT,i,1);
      norm12 = g.nuE0T(nT,i,2);
      norm21 = 0;
      norm22 = 0;
      for j = 1:3
        for nT2 = 1:g.numT
          if markE0TL(nT2,j) || markE0TOS(nT2,j)
            if g.V0E(g.E0T(nT2,j),2) == p1
              norm21 = g.nuE0T(nT2,j,1);
              norm22 = g.nuE0T(nT2,j,2);
              infDir(p1,1) = 0.5 * (norm11 + norm21);
              infDir(p1,2) = 0.5 * (norm12 + norm22);
              nrm = sqrt(infDir(p1,1).^2 + infDir(p1,2).^2);
              infDir(p1,1) = infDir(p1,1)/nrm;
              infDir(p1,2) = infDir(p1,2)/nrm;
              break;
            elseif g.V0E(g.E0T(nT2,j),1) == p2
              norm21 = g.nuE0T(nT2,j,1);
              norm22 = g.nuE0T(nT2,j,2);
              infDir(p2,1) = 0.5 * (norm11 + norm21);
              infDir(p2,2) = 0.5 * (norm12 + norm22);
              nrm = sqrt(infDir(p2,1).^2 + infDir(p2,2).^2);
              infDir(p2,1) = infDir(p2,1)/nrm;
              infDir(p2,2) = infDir(p2,2)/nrm;
            end
          end
        end
      end
    end
  end
end
end%function