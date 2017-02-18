%
function ret = integrateRefEdgeMuPhiInt(N, Nhybrid, basesOnQuad, basesOnGamma)
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = Nhybrid-1;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

ret = zeros(Nhybrid, N, 3, 2); % [N x N]
for nn = 1 : 3
    for i = 1 : Nhybrid
        for j = 1 : N
            ret(i, j, nn, 1) =  W * (basesOnQuad.phi1D{qOrd}(:,j, nn) .* basesOnGamma.thetaPhi1D{qOrd}(:,i, 1) );
            ret(i, j, nn, 2) =  W * (basesOnQuad.phi1D{qOrd}(:,j, nn) .* basesOnGamma.thetaPhi1D{qOrd}(:,i, 2) );
        end % for
    end % for
end % for
end % function
