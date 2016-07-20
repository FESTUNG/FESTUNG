%% assembles global open sea boundary edge contributions of \int_E \varphi_i * ( u * nu^{1} + v * nu^{2} )
%% Dirichlet!


function globUOS = assembleGlobUOSInf(gInf, u, v, HOSAlg, NInf, p, beta, basesOnQuad)
% global gPhiInf1D
KInf = gInf.numInfElem;
pInf = NInf/(p+1)-1;  qOrdInf = max(2*pInf+1, 1); [~, WInf] = quadRuleInf1D(qOrdInf, beta);
qOrd = max(2*p+1, 1); [Q, W] = quadRule1D(qOrd);  %R = length(W);
wDiff = length(WInf)-length(W);
W = [W zeros(1,wDiff)];


Q2X1 = @(X1, X2) gInf.xHat(:,1)*ones(size(X1)) + cos(gInf.Alpha)*X1 - ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (sin(gInf.Beta)*(X2-0.5));
Q2X2 = @(X1, X2) gInf.xHat(:,2)*ones(size(X2)) + sin(gInf.Alpha)*X1 + ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (cos(gInf.Beta)*(X2-0.5));

globUOS = cell(2,1); globUOS{1} = zeros(KInf, NInf);  globUOS{2} = zeros(KInf, NInf);
for n = 2 : 2 %2 entspricht endlicher Kante
  [Q1, Q2] = gammaMapInf(n, Q);
  uDn = u(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
  uDn = [uDn zeros(size(uDn,1), wDiff)];
  
  vDn = v(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
  vDn = [vDn zeros(size(vDn,1), wDiff)];
  
  HDn = HOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
  HDn = [HDn zeros(size(HDn,1), wDiff)];
  Jkn = gInf.detE0I(:,n);
  for i = 1 : NInf
    uu = uDn.*uDn.*HDn * ( W .* basesOnQuad.gPhiInf1D{qOrdInf}(:, i, n)' )';
    uv = uDn.*vDn.*HDn * ( W .* basesOnQuad.gPhiInf1D{qOrdInf}(:, i, n)' )';
    vv = vDn.*vDn.*HDn * ( W .* basesOnQuad.gPhiInf1D{qOrdInf}(:, i, n)' )';
    globUOS{1}(:,i) = globUOS{1}(:,i) + Jkn .* gInf.nuE0Inf(:,n,1) .* uu;
    globUOS{1}(:,i) = globUOS{1}(:,i) + Jkn .* gInf.nuE0Inf(:,n,2) .* uv;  
    
    globUOS{2}(:,i) = globUOS{2}(:,i) + Jkn .* gInf.nuE0Inf(:,n,1) .* uv;
    globUOS{2}(:,i) = globUOS{2}(:,i) + Jkn .* gInf.nuE0Inf(:,n,2) .* vv;
  end % for
end % for
globUOS{1} = reshape(globUOS{1}',KInf*NInf,1);  globUOS{2} = reshape(globUOS{2}',KInf*NInf,1);
end % function
