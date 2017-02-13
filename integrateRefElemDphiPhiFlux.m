%%TODO This needs a better name
%
% I precompute matrices G_bar that allows for an 'easy' evaluation of 
% u_{m} phi_{kj} \partial_{x_m} phi_{ki}
%
function ret = integrateRefElemDphiPhiFlux(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [~,~,W] = quadRule2D(qOrd);
Nip = size(W,2);
ret = zeros( Nip, N, N, 2 ); % [ N x N x 2]

% ret = cell(2,1); 
% ret{1} = sparse(N, N, Nip); 
% ret{2} = sparse(N, N, Nip);

% if N > 1 % p > 0
    for i = 1 : N
        for j = 1 : N
            for m = 1 : 2
                for ip = 1:Nip
%                     W(ip)
%                     basesOnQuad.gradPhi2D{qOrd}(ip,i,m)
%                     basesOnQuad.phi2D{qOrd}(ip,j)
%                     basesOnQuad.gradPhi2D{qOrd}(ip,i,m) .* basesOnQuad.phi2D{qOrd}(ip,j)
                    ret( ip, i, j, m ) =  W(ip) .* basesOnQuad.gradPhi2D{qOrd}(ip,i,m) .* basesOnQuad.phi2D{qOrd}(ip,j) ;
                    ret( ip, i, j, m )
                end
            end % for
        end % for
    end % for
% end % if
end % function
