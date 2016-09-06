function basesOnQuad = computeBasesOnQuad( p , ord )

numBases = (p+1)^2;

[Q1 , Q2 , ~] = gaussQuadRule2D(ord);
basesOnQuad.gPhi2D = zeros(length(Q1), numBases);
basesOnQuad.gGradPhi2D = zeros(length(Q1), numBases, 2);

[Q, ~] = gaussQuadRule1D(ord);
basesOnQuad.gPhi1D = zeros(length(Q), numBases, 4);

for i = 1 : numBases
    basesOnQuad.gPhi2D(:,i) = phi(i, Q1, Q2);
    basesOnQuad.gGradPhi2D(:,i,1) = gradPhi(i, 1, Q1, Q2);
    basesOnQuad.gGradPhi2D(:,i,2) = gradPhi(i, 2, Q1, Q2);
    basesOnQuad.gPhi1D(:,i,1) = phi(i, Q, zeros(length(Q),1));
    basesOnQuad.gPhi1D(:,i,2) = phi(i, Q, ones(length(Q),1));
    basesOnQuad.gPhi1D(:,i,3) = phi(i, ones(length(Q),1), Q);
    basesOnQuad.gPhi1D(:,i,4) = phi(i, zeros(length(Q),1), Q);
end  % for

basesOnQuad.gPsi1D = zeros(length(Q), p);
basesOnQuad.gDerivPsi1D = zeros(length(Q), p);
basesOnQuad.gPsi1Dbnd = zeros(2, p);

for i = 1 : p+1
    basesOnQuad.gPsi1D(:,i) = psi1D(i, Q);
    basesOnQuad.gDerivPsi1D(:,i) = derivPsi(i, Q);
    basesOnQuad.gPsi1Dbnd(1,i) = psi1D(i,1);
    basesOnQuad.gPsi1Dbnd(2,i) = psi1D(i,0);
end  % for

end  % function