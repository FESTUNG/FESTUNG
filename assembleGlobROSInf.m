function globROS = assembleGlobROSInf(gInf, HOSAlg, NInf, p, gConst, beta, basesOnQuad)
%% assembles global open sea boundary edge contributions of -0.5 * gConst * \int_E \varphi_i * HOS^2 * nu
%% globROS corresponds to open sea boundary edge integrals in weak form of -0.5 * gConst * HOS^2 in PDE
% global gPhiInf1D
KInf = gInf.numInfElem; pInf = NInf/(p+1)-1;
qOrdInf = 2*pInf+1; qOrd = 2*p+1;  [Q, W] = quadRule1D(qOrd); 
[~, WInf] = quadRuleInf1D(qOrdInf, beta); wDiff = length(WInf)- length(W);
W =[W, zeros(1,wDiff)];
Q2X1 = @(X1, X2) gInf.xHat(:,1)*ones(size(X1)) + cos(gInf.Alpha)*X1 - ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (sin(gInf.Beta)*(X2-0.5));
Q2X2 = @(X1, X2) gInf.xHat(:,2)*ones(size(X2)) + sin(gInf.Alpha)*X1 + ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (cos(gInf.Beta)*(X2-0.5));
globROS = cell(2, 1);
globROS{1} = zeros(KInf, NInf); globROS{2} = zeros(KInf, NInf);
for n = 2 : 2
  [Q1, Q2] = gammaMapInf(n, Q);
  HOSn = HOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
  HOSn =[HOSn, zeros(size(HOSn,1),wDiff)];
  HOSkn = -gInf.detE0I(:,n);
  HOSnWphi = bsxfun(@times, HOSn .* HOSn, W) * basesOnQuad.gPhiInf1D{qOrdInf}(:, :, n);
  globROS{1} = globROS{1} + bsxfun(@times, HOSnWphi, HOSkn .* gInf.nuE0Inf(:,n,1));
  globROS{2} = globROS{2} + bsxfun(@times, HOSnWphi, HOSkn .* gInf.nuE0Inf(:,n,2));
end % for
globROS{1} = reshape(0.25 * gConst * globROS{1}.',KInf*NInf,1);  %Average: 0.5 -> 0.25
globROS{2} = reshape(0.25 * gConst * globROS{2}.',KInf*NInf,1);
end % function
