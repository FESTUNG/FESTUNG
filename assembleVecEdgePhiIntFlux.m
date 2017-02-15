function ret = assembleVecEdgePhiIntFlux( g, N, fluxEval, markE0Tbdr, basesOnQuad )
%Assert
K = g.numT;

validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

% KEdge = g.numE;

% Assemble matrix
ret = zeros( K*N, 1 );
% warning('flipping in assembleVecEdgePhiIntFlux');

myI = 1;

for iT = 1:K
    for iE = 1:3
        edgeNr = g.E0T(iT, iE);
        flip = 1;
        
        fluxLocal = fluxEval( :, :, iT, iE);
        
        bases = basesOnQuad.phi1D{qOrd}( :, :, iE);
        if (g.T0E(edgeNr, 2) == iT)
            flip = 2;
            %             bases = flipud(bases);
%             fluxLocal = flipud( fluxLocal );
        end
%         if (2 == iT)
% %             flip = 2;
% %             bases = flipud(bases);
% %              fluxLocal = flipud( fluxLocal );
%         end
        
        
        iTs = (iT-1)*N + 1;
        iTe = (iT)*N;
        
        tmp = zeros(N, 1);
        fluxEvalNu = zeros( 1, size(W,2));
        for m=1:2
            fluxEvalNu =  fluxEvalNu + fluxLocal( m, : ) .* g.nuE0T( iT, iE, m );
            %             fluxEvalNu =  fluxEvalNu + fluxEval( m, :, iT, iE ) .* g.nuE0T( iT, iE, m );
            %             fluxEvalNu =  fluxEvalNu + fliplr( fluxEval( m, :, edgeNr ) ) .* g.nuE0T( iT, iE, m );
            
        end
        
        tmpAsdf( :, iT, iE ) = fluxEvalNu;
        for i = 1:N
            %             g.areaE( edgeNr )
            %             markE0Tbdr(iT, iE)
            %             fluxEvalNu'
            %             basesOnQuad.phi1D{qOrd}( :, i, iE)
            %             ( W * (fluxEvalNu' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) )
            %             g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * (fluxEvalNu' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) )
            %             tmp(i) = tmp(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * (fluxEvalNu' .* basesOnQuad.phi1D{qOrd}( :, i, iE) ) );
%             g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * (fluxEvalNu' .* bases( :, i ) ) )
            tmp(i) = tmp(i) + g.areaE( edgeNr ) .* markE0Tbdr(iT, iE) .* ( W * (fluxEvalNu' .* bases( :, i ) ) );
        end

        qwert( myI:myI+N-1 ) = tmp(:);
        myI = myI + N;
        
        
        ret(iTs:iTe) = ret(iTs:iTe) + tmp;
    end
end


end