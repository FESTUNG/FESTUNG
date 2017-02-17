function ret = assembleVecEdgePhiIntVal( g, N, cEval, markE0Tbdr, problemData, basesOnQuad )
%Assert
K = g.numT;
Kedge = g.numE;

validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

% warning('flipping in assembleVecEdgePhiIntVal');
% warning('additional return for testing in assembleVecEdgePhiIntVal');

% warning('Setting Dirichlet boundary on every edge assembleVecEdgePhiIntVal');


% Assemble matrix
ret = zeros( K*N, 1 );
% ret2 = zeros( K*N, 1 );
% for iT = 1:K
%     for iE = 1:3
%         edgeNr = g.E0T(iT, iE);
%         iTs = (iT-1)*N + 1;
%         iTe = (iT)*N;
%
%         tmp = zeros(N, 1);
%         %         tmp2 = zeros(N, 1);
%         if ( g.markE0TbdrD(iT, iE) )
%
%             for i = 1:N
%                 %             g.areaE( edgeNr )
%                 %             markE0Tbdr(iT, iE)
%                 %             W
%                 %             cEval(edgeNr,:)'
%                 %             basesOnQuad.phi1D{qOrd}( :, i, iE)
%                 %             ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) )
%                 %             ( W * ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) )
%                 %             tmp(i) = tmp(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * ( fliplr(cEval(edgeNr,:))' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%
%
%                 %Set Dirichlet on every edge
%                 %             tmp(i) = tmp(i) + g.areaE( edgeNr ) .* ( W * ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%
%                 %Actual implementation
% %                 tmp(i) = tmp(i) + g.areaE( edgeNr ) .* g.markE0TbdrD(iT, iE) .* ( W * ( cEval(edgeNr,:)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%                 tmp(i) = tmp(i) + g.areaE( edgeNr ) .* g.markE0TbdrD(iT, iE) .* ( W * ( cEval(iT,:, iE)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%                 %             tmp2(i) = tmp2(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * ( fliplr(cEval(edgeNr,:))' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%
%             end
%             %         if (iT == 2 && edgeNr==3)
%             %             ret(iTs:iTe) = ret(iTs:iTe) + fliplr( tmp );
%             %         else
%             %             ret(iTs:iTe) = ret(iTs:iTe) + tmp;
%             %         end
% %         else
% %             bases = problemData.basesOnGamma.phi1D{qOrd};
% %             cDiscLambda = problemData.cDiscLambda;
% %
% %             lambdaLocal = cDiscLambda(edgeNr,:) * bases';
% %             for i = 1:N
% %                 %Actual implementation
% %                 tmp(i) = tmp(i) + g.areaE( edgeNr ) .* g.markE0TbdrN(iT, iE) .* ( W * ( lambdaLocal(:) .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
% %                 %             tmp2(i) = tmp2(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * ( fliplr(cEval(edgeNr,:))' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
% %
% %             end
%         end
%
%         %Actual implementation
%         ret(iTs:iTe) = ret(iTs:iTe) + tmp;
%
%         %         ret2(iTs:iTe) = ret2(iTs:iTe) + tmp2;
%     end
% end

for n=1:3
%     tmp(i) = tmp(i) + g.areaE( edgeNr ) .* g.markE0TbdrD(iT, iE) .* ( W * ( cEval(iT,:, iE)' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
    fac = g.areaE0T( :, n ) .* g.markE0TbdrD(:, n);
%     ret(:) = ret(:) +  fac .* ( W * ( cEval(:,:, n)' .* basesOnQuad.phi1D{qOrd}( :, :, n) ) );
    ret(:) =ret(:) + reshape( ( fac .* cEval(:,:, n) * (W' .* basesOnQuad.phi1D{qOrd}( :, :, n) ))', K*N, 1 );
end

end
