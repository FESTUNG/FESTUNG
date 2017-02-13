function ret = assembleVecEdgePhiIntVal( g, N, cEval, markE0Tbdr, basesOnQuad )
%Assert
K = g.numT;
Kedge = g.numE;

validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

% warning('flipping in assembleVecEdgePhiIntVal');
% warning('additional return for testing in assembleVecEdgePhiIntVal');

% Assemble matrix
ret = zeros( K*N, 1 );
% ret2 = zeros( K*N, 1 );
for iT = 1:K
    for iE = 1:3
        edgeNr = g.E0T(iT, iE);
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        
        tmp = zeros(N, 1);
%         tmp2 = zeros(N, 1);

        for i = 1:N
%             g.areaE( edgeNr )
%             markE0Tbdr(iT, iE)
%             W
%             cEval(edgeNr,:)'
%             basesOnQuad.phi1D{qOrd}( :, i, iE) 
%             ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) )
%             ( W * ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) )
%             tmp(i) = tmp(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * ( fliplr(cEval(edgeNr,:))' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
            tmp(i) = tmp(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%             tmp2(i) = tmp2(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * ( fliplr(cEval(edgeNr,:))' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
        end
        
        ret(iTs:iTe) = ret(iTs:iTe) + tmp;
%         ret2(iTs:iTe) = ret2(iTs:iTe) + tmp2;
    end
end
end