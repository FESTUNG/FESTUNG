function ret = assembleVecV0T1DPhiIntFuncDiscIntNuHeight(g1D, markV0T, dataDisc, heightV0T1D, barN, qOrd, basesOnQuad)
barK = g1D.numT;
ret = zeros(barK, barN);
for n = 1 : 2
  funcDiscV0T = (dataDisc{1} + (n-1) * dataDisc{2}) * basesOnQuad.phi0D{qOrd}(:, n);
  ret = ret + (markV0T(:, n) .* g1D.nuV0T(:, n) .* funcDiscV0T ./ heightV0T1D(:, n)) * basesOnQuad.phi0D{qOrd}(:, n).';
end  % for n
ret = reshape(ret.', barK*barN, 1);
end  % function