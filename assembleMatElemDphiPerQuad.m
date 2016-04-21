function ret = assembleMatElemDphiPerQuad(g, refElemDphiPerQuad)
K = g.numT;
ret = cell(2,1);
ret{1} = kron(spdiags(g.B(:,2,2), 0, K, K), refElemDphiPerQuad(:,:,1)) - kron(spdiags(g.B(:,2,1), 0, K, K), refElemDphiPerQuad(:,:,2));
ret{2} = kron(spdiags(g.B(:,1,1), 0, K, K), refElemDphiPerQuad(:,:,2)) - kron(spdiags(g.B(:,1,2), 0, K, K), refElemDphiPerQuad(:,:,1));
end % function
