%
function ret = integrateRefEdgeMuPhiIntFlux(N, Nlambda, basesOnQuad, basesOnGamma)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);
Nip = size(W,2);
%ret = zeros(Nip, N, Nlambda, 3, 2); % [N x N x 3]
ret = cell(2,1);
ret{1} = zeros(N, Nlambda, 3, Nip); % [N x N x 3]
ret{2} = zeros(N, Nlambda, 3, Nip); % [N x N x 3]
for n = 1 : 3 % 3 edges
    for i = 1 : N
        for j = 1 : Nlambda
            for ip = 1:Nip
                ret{1}(i, j, n, ip) = W(ip) .* basesOnQuad.phi1D{qOrd}(ip,i,n) .* basesOnGamma.phi1D{qOrd}(ip,j) ;
                baseFlip = flipud(basesOnGamma.phi1D{qOrd});
                ret{2}(i, j, n, ip) = W(ip) * (basesOnQuad.phi1D{qOrd}(ip,i,n) .* baseFlip(ip,j) );
            end
        end % for
    end % for
end % for
end % function
