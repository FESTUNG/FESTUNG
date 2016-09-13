function globM = assembleMatElemTrapPhiPhi(g, refElemPhiPhi)
globM = kron(spdiags(g.detJ0T{1}, 0, g.numT, g.numT), refElemPhiPhi{1}) + ...
        kron(spdiags(g.detJ0T{2}, 0, g.numT, g.numT), refElemPhiPhi{2});
end  % function