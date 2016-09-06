function computeBasesOnQuad( p , ord )

global gPhi2D gGradPhi2D gPhi1D gPsi1D gDerivPsi1D gPsi1Dbnd

numBases = (p+1)^2;

[Q1 , Q2 , ~] = gaussQuadRule2D(ord);
gPhi2D = zeros(length(Q1), numBases);
gGradPhi2D = zeros(length(Q1), numBases, 2);

[Q, ~] = gaussQuadRule1D(ord);
gPhi1D = zeros(length(Q), numBases, 4);

for i = 1 : numBases
    gPhi2D(:,i) = phi(i, Q1, Q2);
    gGradPhi2D(:,i,1) = gradPhi(i, 1, Q1, Q2);
    gGradPhi2D(:,i,2) = gradPhi(i, 2, Q1, Q2);
    gPhi1D(:,i,1) = phi(i, Q, zeros(length(Q),1));
    gPhi1D(:,i,2) = phi(i, Q, ones(length(Q),1));
    gPhi1D(:,i,3) = phi(i, ones(length(Q),1), Q);
    gPhi1D(:,i,4) = phi(i, zeros(length(Q),1), Q);
end  % for

gPsi1D = zeros(length(Q), p);
gDerivPsi1D = zeros(length(Q), p);
gPsi1Dbnd = zeros(2, p);

for i = 1 : p+1
    gPsi1D(:,i) = psi1D(i, Q);
    gDerivPsi1D(:,i) = derivPsi(i, Q);
    gPsi1Dbnd(1,i) = psi1D(i,1);
    gPsi1Dbnd(2,i) = psi1D(i,0);
end  % for


end  % function