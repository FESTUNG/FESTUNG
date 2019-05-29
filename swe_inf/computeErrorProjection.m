function [errH, erruH, errvH] = computeErrorProjection(gInf, HAlg, uAlg, vAlg, cDGInf, p, pInf, beta)
%qOrd = max(2*p, 1); qOrdInf = max(2*pInf,1); % qOrdInf = 126;
qOrd =11; qOrdInf=32;
NInf = (p+1)*(pInf+1);
[Q1, Q2, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
R = length(W);

uHAlg = @(x,y) uAlg(x,y).*HAlg(x,y); 
vHAlg = @(x,y) vAlg(x,y).*HAlg(x,y); 



F1 = @(X1, X2) gInf.xHat(:,1)*ones(size(X1)) + cos(gInf.Alpha)*X1 - ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (sin(gInf.Beta)*(X2-0.5));
F2 = @(X1, X2) gInf.xHat(:,2)*ones(size(X2)) + sin(gInf.Alpha)*X1 + ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (cos(gInf.Beta)*(X2-0.5));

detQuadP = zeros(gInf.numInfElem,size(Q1,2));
for nq = 1 : size(Q1,2)
  detQuadP(:,nq) = (cos(gInf.Alpha)-sin(gInf.Beta).*2.*gInf.M.*(Q2(nq)-0.5)).*cos(gInf.Beta).*(gInf.H + 2*Q1(nq)*gInf.M) + ...
    sin(gInf.Beta).*(gInf.H+2*Q1(nq)*gInf.M).*(sin(gInf.Alpha) + 2*cos(gInf.Beta).*gInf.M*(Q2(nq)-0.5));
end%for

for i = 1:3
  err=0;
  switch i
    case 1, fAlg = @(x,y) HAlg(x,y);
    case 2, fAlg = @(x,y) uHAlg(x,y);
    case 3, fAlg = @(x,y) vHAlg(x,y);
  end % switch
  for r = 1:R
    cDGOnQuadPts = 0;
    for n = 1:NInf
      cDGOnQuadPts = cDGOnQuadPts + cDGInf(:,n,i) * phiInf(n, p, Q1(r), Q2(r), beta);
    end%for
    fExOnQuadPts = fAlg(F1(Q1(r),Q2(r)),F2(Q1(r),Q2(r)));
    err = err + W(r).*(fExOnQuadPts - cDGOnQuadPts).^2.*exp(-beta*Q1(r)).*detQuadP(:,r);
  end%for
  switch i
    case 1, errH  = sqrt(sum(err));
    case 2, erruH = sqrt(sum(err));
    case 3, errvH = sqrt(sum(err));
  end % switch
end%for
end % function
