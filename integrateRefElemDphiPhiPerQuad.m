%%TODO This needs a better name
%
% I precompute matrices G_bar that allows for an 'easy' evaluation of 
% u_{m} phi_{kj} \partial_{x_m} phi_{ki}
%
function ret = integrateRefElemDphiPhiPerQuad(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [~,~,W] = quadRule2D(qOrd);
R = size(W,2);
ret = cell(2,1);  % {2} x [ N x N x 2]
ret{1} = zeros(N, N, R);
ret{2} = zeros(N, N, R);
for i = 1 : N
    for j = 1 : N
        for m = 1 : 2
            for ip = 1:R
                ret{m}( i, j, ip ) =  W(ip) .* basesOnQuad.phi2D{qOrd}(ip,j) .* basesOnQuad.gradPhi2D{qOrd}(ip,i,m) ;
            end
        end % for
    end % for
end % for
end % function
