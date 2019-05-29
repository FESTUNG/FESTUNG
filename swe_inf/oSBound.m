function markVOSBound = oSBound(g, markE0TOS, markE0TL)
markVOSBound = zeros(g.numV, 1);
for nT = 1: g.numT
  for i = 1:3
    if (markE0TOS(nT, i))
      p1 = g.V0E(g.E0T(nT, i),1);
      p2 = g.V0E(g.E0T(nT, i),2);
    else
      continue;
    end%if
    for nT2 = 1: g.numT
      for j = 1:3
        if (markE0TL(nT2, j))
          q1 = g.V0E(g.E0T(nT2, j),1);
          q2 = g.V0E(g.E0T(nT2, j),2);
          if p2 == q1 
            markVOSBound(p2)=1;
          end%if
          if p1 == q2
            markVOSBound(p1)=1;
          end%if
        end%if
      end%for
    end%for
  end%for
end%for
end%function