%
function ret = integrateRefEdgePhiIntPhiIntFlux(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
Nip = size(W,2);
ret = zeros(N, N, 3, Nip); % [N x N x 3]
for n = 1 : 3 % 3 edges
    for i = 1 : N
        for j = 1 : N
            for ip = 1:Nip
                ret(i, j, n, ip) = W(ip) .* basesOnQuad.phi1D{qOrd}(ip,i,n) .* basesOnQuad.phi1D{qOrd}(ip,j,n) ;
            end
        end % for
    end % for
end % for
end % function
