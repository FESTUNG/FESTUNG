function ret = assembleMatElemPhiPhiFuncDiscConst(g, refElemPhiPhi, dataDiscConst)
K = g.numT;
ret = 2 * kron(spdiags(g.areaT .* dataDiscConst, 0, K, K), refElemPhiPhi);
end % function
