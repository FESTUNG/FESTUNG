%
function ret = integrateRefEdgeMuPhiIntFlux(N, Nmu, basesOnQuad, basesOnGamma)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
Nip = size(W,2);
ret = cell(2,1);
ret{1} = zeros(N, Nmu, 3, Nip); % [N x N x 3]
ret{2} = zeros(N, Nmu, 3, Nip); % [N x N x 3]
for n = 1 : 3 % 3 edges
    for i = 1 : N
        for j = 1 : Nmu
            for ip = 1:Nip
                ret{1}(i, j, n, ip) = W(ip) .* basesOnQuad.phi1D{qOrd}(ip,i,n) ...
                    .* basesOnGamma.thetaPhi1D{qOrd}(ip, j, 1);
                ret{2}(i, j, n, ip) = W(ip) .* basesOnQuad.phi1D{qOrd}(ip,i,n) ...
                    .* basesOnGamma.thetaPhi1D{qOrd}(ip, j, 2);
            end
        end % for
    end % for
end % for
end % function
