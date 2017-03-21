function ret = integrateRefEdgeTrapPhi1DIntPerQuad(N, qOrd, basesOnQuad)
[~, W] = quadRule1D(qOrd); R = length(W);
ret = zeros(N, R, 4);
for n = 1 : 2
  for i = 1 : N
    [m, ~] = execin('darcyVert/mapTensorProductIndexInv', i);
    ret(i,:,n) = W .* basesOnQuad.phi1D{qOrd}(:,m).';
  end  % for i
end  % for n
for n = 3 : 4
  for i = 1 : N
    [m, ~] = execin('darcyVert/mapTensorProductIndexInv', i);
    ret(i,:,n) = W * basesOnQuad.phi0D{qOrd}(m,n-2);
  end  % for i
end  % for n
end  % function