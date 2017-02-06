%
function ret = integrateRefEdgePhiIntMu(N, Nhybrid, basesOnQuad, basesOnGamma)
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = Nhybrid-1;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

ret = zeros(Nhybrid, N, 3); % [N x N]
for nn = 1 : 3
  for i = 1 : Nhybrid
    for j = 1 : N
%         basesOnGamma.phi1D{qOrd}(:,i)
%         basesOnQuad.phi1D{qOrd}(:,j)
%         W'
      ret(i, j, nn) = sum( W' .* basesOnGamma.phi1D{qOrd}(:,i) .* basesOnQuad.phi1D{qOrd}(:,j, nn) );
    end % for
  end % for
end % for  
end % function
