function plotFunction ( cDG, gInf, beta)
global gPhiInf2D
[~, NInf] = size(cDG(:,:,1));
p = sqrt(NInf)-1;
qOrd = max(2*p,1);
[Q1, Q2, ~] = quadRuleInf2D(qOrd, beta, 2);
R = length(Q1);


X1 = kron(gInf.xHat(:,1),ones(1,R))+kron(cos(gInf.Alpha),Q1)+kron(-sin(gInf.Beta).*gInf.H,Q2-0.5)+kron(-sin(gInf.Beta)*2.*gInf.M, Q1.*Q2-0.5*Q1);
X2 = kron(gInf.xHat(:,1),ones(1,R))+kron(sin(gInf.Alpha),Q1)+kron( cos(gInf.Beta).*gInf.H,Q2-0.5)+kron( cos(gInf.Beta)*2.*gInf.M, Q1.*Q2-0.5*Q1);

hApprxOnQuadPts = cDG(:,:,1)*gPhiInf2D{qOrd}(:,:,2)';
uApprxOnQuadPts = cDG(:,:,2)*gPhiInf2D{qOrd}(:,:,2)';
vApprxOnQuadPts = cDG(:,:,3)*gPhiInf2D{qOrd}(:,:,2)';

size(X1)
size(X2)
size(hApprxOnQuadPts)
%[x,y,z] = meshgrid(X1,X2,hApprxOnQuadPts);

surf(X1,X2, hApprxOnQuadPts);
surf(X1,X2, uApprxOnQuadPts);
surf(X1,X2, vApprxOnQuadPts);
end % function