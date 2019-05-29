function globVOSRiemInf = assembleGlobVOSRiemInf(gInf, HOSAlg, uOSAlg, vOSAlg, PhiPhi1D, markE0TOS, cDGInf, gConst, p, beta, gCdiagInf, gCoffdiagInf, basesOnQuad)
%% assembles global open sea boundary edge contributions of \int_E \varphi_i * \lambda * \varphi_j
%% by building a matrix and 
%% assembles global open sea boundary edge contributions of \int_E \varphi_i * \lambda * HOS
%% by building a vector
%% here \lambda = abs( ( (H^0.5 * uH)^- + (H^0.5 * uH)^+ ) * nu^{1} + ( (H^0.5 * vH)^- + (H^0.5 * vH)^+ ) * nu^{2} ) / ((H^-)^1.5 + (H^+)^1.5)
%%              + ( 0.5 * gConst * ((H^-)^1.5 + (H^+)^1.5) / ((H^-)^0.5 + (H^+)^0.5) )^0.5
%% with H^- being the unknown height and H^+ the perscribed open sea boundary height
%% and uH^- = uH^+, vH^- = vH^+ respectively being the values of the unknown momenta from the interior<-------changed!
%% multiplying globVOSRiem{1} with H corresponds to open sea boundary edge integrals resulting from using a Lax-Friedrichs solver
%% globVOSRiem{2} corresponds to open sea boundary edge integrals resulting from using a Lax-Friedrichs solver

%%---->geaendert!
% global gPhiInf1D gCdiagInf gCoffdiagInf
[KInf, NInf, ~] = size(cDGInf);
globVOSRiemInf = cell(2,1);
globVOSRiemInf{1} = sparse(KInf*NInf, KInf*NInf); 
globVOSRiemInf{2} = zeros(KInf,NInf);
pInf = NInf/(p+1)-1;  qOrdInf = max(2*pInf+1, 1); [~, WInf] = quadRuleInf1D(qOrdInf, beta);
qOrd = max(2*p+1, 1); [Q, W] = quadRule1D(qOrd); R = length(W);
Q2X1 = @(X1, X2) gInf.xHat(:,1)*ones(size(X1)) + cos(gInf.Alpha)*X1 - ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (sin(gInf.Beta)*(X2-0.5));
Q2X2 = @(X1, X2) gInf.xHat(:,2)*ones(size(X2)) + sin(gInf.Alpha)*X1 + ...
  (gInf.H*ones(size(X1)) + 2*gInf.M*X1).* (cos(gInf.Beta)*(X2-0.5));
lambda = zeros(KInf,length(WInf));

wDiff = length(WInf)-R;
W = [W, zeros(1,wDiff)];
[Q1, Q2] = gammaMapInf(2, Q);
HOS = HOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
uOS = uOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
vOS = vOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
HOS = [HOS, zeros(size(HOS,1),wDiff)];
uOS = [uOS, zeros(size(HOS,1),wDiff)];
vOS = [vOS, zeros(size(HOS,1),wDiff)];
for n = 1 : 3
  for r = 1 : R
    HLPow1o2 = gCdiagInf{1,2}(:,r).^0.5;
    HRPow1o2 = HOS(:,r).^0.5;
    lambda(:,r) =   abs( ( ( HRPow1o2 .* gCdiagInf{2,2}(:,r) + HLPow1o2 .* uOS(:,r) ) .* gInf.nuE0Inf(:,2,1)   ...
                         + ( HRPow1o2 .* gCdiagInf{3,2}(:,r) + HLPow1o2 .* vOS(:,r) ) .* gInf.nuE0Inf(:,2,2) ) ...
                         ./( HRPow1o2 .* gCdiagInf{1,2}(:,r) + HLPow1o2 .* HOS(:,r)) );
    lambda(:,r) = ( min(lambda(:,r), 0) + max(lambda(:,r), 0) ...
             + sqrt( gConst * (gCdiagInf{1,2}(:,r).^1.5 + HOS(:,r).^1.5) ./ (HLPow1o2 + HRPow1o2) ) );
    if n == 2
        globVOSRiemInf{1} = globVOSRiemInf{1} + 0.5 * kron( spdiags( gInf.detE0I(:,n) .* lambda(:,r), 0, KInf, KInf ), PhiPhi1D(:,:,r,n) );
    end % if
  end % for
  globVOSRiemInf{2} = globVOSRiemInf{2} + ( repmat(markE0TOS(gInf.T0I,n),1,length(WInf)) .* HOS .* lambda .* ( gInf.detE0I(:,2) * W ) ) * basesOnQuad.gPhiInf1D{qOrdInf}(:, :, 2);
end % for
globVOSRiemInf{2} = reshape(0.5 * globVOSRiemInf{2}.',KInf*NInf,1);
end % function

% for n = 2 : 2
%   [Q1, Q2] = gammaMapInf(n, Q);
%   HOSn = HOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
%   uOSn = uOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
%   vOSn = vOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2));
%   for r = 1:R
%     HPow1o2l =  gCdiagInf{1,n}(:,r).^0.5;
%     HPow1o2r =  HOSn(:,r).^0.5;
%     lambda(:,r) =   abs( ( ( HPow1o2r.*gCdiagInf{2,n}(:,r) + HPow1o2l.*uOSn(:,r).*HOSn(:,r) ) .* gInf.nuE0Inf(:,n,1)    ...
% 												 + ( HPow1o2r.*gCdiagInf{3,n}(:,r) + HPow1o2l.*vOSn(:,r).*HOSn(:,r) ) .* gInf.nuE0Inf(:,n,2) )  ...
% 											   ./( gCdiagInf{1,n}(:,r).*HPow1o2r + HOSn(:,r).*HPow1o2l ) );
%     lambda(:,r) = min(lambda(:,r), 0) + max(lambda(:,r), 0); ...
%            + sqrt( gConst * (gCdiagInf{1,n}(:,r).^1.5 + HOSn(:,r).^1.5) ./ (HPow1o2l + HPow1o2r) );
%     globVOSRiem{1} = globVOSRiem{1} + 0.5 * kron( spdiags( gInf.detE0I(:,n) .* lambda(:,r), 0, KInf, KInf ), PhiPhi1D(:,:,r,n) );
%   end % for
%   globVOSRiem{2} = globVOSRiem{2} + ( HOSn .* lambda .* ( gInf.detE0I(:,n) * W ) ) * gPhiInf1D{qOrd}(:, :, n);
% end % for
