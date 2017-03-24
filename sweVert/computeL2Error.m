function [ error ] = computeL2ErrorSub( g , fDG , fAlg , ord )

global gPhi2D

[Q1, Q2, W] = gaussQuadRule2D(ord);
QX = kron(g.deltaX*ones(g.numTsub,1),Q1') + kron(g.coordV0Tsub(:,1,1),ones(size(Q1')));
QY = kron(g.BAsub, Q1') + kron(g.DAsub, Q2') + kron(g.ACBDsub, (Q1 .* Q2)')...
    + kron(g.coordV0Tsub(:,1,2), ones(size(Q1')));

error = sqrt( dot( (fDG * gPhi2D' - fAlg(QX, QY)).^2  * W' , g.areaTsub ) );

end  % function