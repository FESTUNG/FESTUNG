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

% warning('gebastel in assembleMatEdgePhiIntPhiIntHybrid. Ich schalte das mal ab, weil ueberalle Dirichlet BC gelten');
% warning('gebastel in assembleMatEdgePhiIntPhiIntHybrid. Ich schalte es ueberall ein');

ret = sparse(K*N, K*N);
% for iT = 1:K
%     iTs = (iT-1)*N + 1;
%     iTe = (iT)*N;
%     for iE = 1:3
%         % Last used version
%         ret( iTs:iTe,  iTs:iTe  ) = ret( iTs:iTe,  iTs:iTe ) + g.markE0Tint(iT, iE) .* g.areaE( g.E0T(iT, iE) ) .* refEdgePhiIntPhiInt(:,:,iE);
%     end
% end

%TODO: ALEX THIS HERE IS THE OPTIMIZED VERSION
% ret2 = sparse(K*N, K*N);
for n = 1 : 3
    % ret2 = ret2 + kron(spdiags( g.areaE( g.E0T(:, n) ),0,K,K ), refEdgePhiIntPhiInt(:,:,n))  ;
    ret( :, : ) = ret( :, : ) + kron(spdiags( g.markE0Tint(:, n) .* g.areaE( g.E0T(:, n) ) ,0,K,K), refEdgePhiIntPhiInt(:,:,n) );
end % for

end % function
