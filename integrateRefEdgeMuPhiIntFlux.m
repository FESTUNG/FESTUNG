%
function ret = integrateRefEdgeMuPhiIntFlux(N, Nlambda, basesOnQuad, basesOnGamma)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
Nip = size(W,2);
ret = zeros(Nip, N, Nlambda, 3, 2); % [N x N x 3]
for n = 1 : 3 % 3 edges
    for i = 1 : N
        for j = 1 : Nlambda
            for ip = 1:Nip
                ret(ip, i, j, n, 1) = W(ip) .* basesOnQuad.phi1D{qOrd}(ip,i,n) .* basesOnGamma.phi1D{qOrd}(ip,j) ;
                baseFlip = flipud(basesOnGamma.phi1D{qOrd});
                ret(ip, i, j, n, 2) = W(ip) * (basesOnQuad.phi1D{qOrd}(ip,i,n) .* baseFlip(ip,j) );
            end
        end % for
    end % for
end % for
end % function
