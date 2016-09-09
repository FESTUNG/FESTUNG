function globM = assembleMatElemTrapPhiPhi(g, hatM)
globM = kron(spdiags(g.deltaJ0T{1}, 0, g.numT, g.numT), hatM{1}) + ...
        kron(spdiags(g.deltaJ0T{2}, 0, g.numT, g.numT), hatM{2});
end  % function