% TODO Check this
% Compute global mas matrix of lambda
% If it is an internal edge, the entry has to be scaled by 2*alpha
% check the element transformation
function ret = assembleMatEdgeOfElemMuMu(g, markE0T, refEdgeMuMu)

K = g.numT;  Nlambda = size(refEdgeMuMu, 1);
KEdge = g.numE;
ret = sparse(KEdge*Nlambda, KEdge*Nlambda);

%Interior edges
ret = sparse(KEdge*Nlambda, KEdge*Nlambda);
for n = 1:3
    Kkn = g.areaE0T( :, n ) .*  markE0T(:, n) ;
    ret = ret + kron( sparse( g.E0T(:, n), g.E0T(:, n), Kkn, KEdge, KEdge ), refEdgeMuMu );
end

% %Exterior edges
% for n = 1:3
%     Kkn = g.areaE0T( :, n ) .*  ( ~g.markE0Tint(:, n) ) ;
%     ret = ret + kron( sparse( g.E0T(:, n), g.E0T(:, n), Kkn, KEdge, KEdge ), refEdgeMuMu );
% end
% 
% %Neumann boundary
% for n = 1:3
%     Kkn = g.areaE0T( :, n ) .*   g.markE0TbdrN(:, n) ;
%     ret = ret - kron( sparse( g.E0T(:, n), g.E0T(:, n), Kkn, KEdge, KEdge ), refEdgeMuMu );
% end
end % function
