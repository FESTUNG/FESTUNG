%
function ret = assembleMatEdgePhiIntPhiIntHybrid(g, refEdgePhiIntPhiInt)
% Extract dimensions
K = g.numT;  N = size(refEdgePhiIntPhiInt, 1);

% Check function arguments that are directly used
validateattributes(refEdgePhiIntPhiInt, {'numeric'}, {'size', [N N 3]}, mfilename, 'refEdgePhiIntPhiInt');

% Assemble matrix
% ret = sparse(K*N, K*N);
% for n = 1 : 3
%   ret = ret + kron(spdiags(g.areaE,0,K,K), refEdgePhiIntPhiInt(:,:,n));
% end % for

ret = sparse(K*N, K*N);
for iT = 1:K
    for iE = 1:3
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        iEs = (g.E0T(iT, iE) - 1)*N + 1;
        iEe = (g.E0T(iT, iE))*N;
        ret( iTs:iTe,  iTs:iTe  ) = ret( iTs:iTe,  iTs:iTe ) + g.areaE( g.E0T(iT, iE) ) .* refEdgePhiIntPhiInt(:,:,iE);
%         ret = ret + kron(spdiags(g.areaE,0,K,K), refEdgePhiIntPhiInt(:,:,n));
    end
end

%TODO: ALEX THIS HERE IS THE OPTIMIZED VERSION
% ret2 = sparse(K*N, K*N);
% for n = 1 : 3
% ret2 = ret2 + kron(spdiags( g.areaE( g.E0T(:, n) ),0,K,K ), refEdgePhiIntPhiInt(:,:,n))  ;
% end % for

end % function
