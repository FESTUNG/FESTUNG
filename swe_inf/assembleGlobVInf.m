function globV = assembleGlobVInf(gInf, PhiPhidiag, PhiPhioffdiag, cDG, gConst, p, beta, gCdiagInf, gCoffdiagInf)
%% assembles global interior edge contributions of \int_E \varphi_i * \lambda * \varphi_j
%% here \lambda = abs( ( (H^0.5 * uH)^- + (H^0.5 * uH)^+ ) * nu^{1} + ( (H^0.5 * vH)^- + (H^0.5 * vH)^+ ) * nu^{2} ) / ((H^-)^1.5 + (H^+)^1.5)
%%              + ( 0.5 * gConst * ((H^-)^1.5 + (H^+)^1.5) / ((H^-)^0.5 + (H^+)^0.5) )^0.5
%% multiplying globV with H (in the continuity equation), uH (in the first momentum equation), vH (in the second momentum equation)
%% corresponds to interior edge integrals resulting from using a Lax-Friedrichs solver
% global gCdiagInf gCoffdiagInf
[KInf, NInf, ~] = size(cDG);
globV = sparse(KInf*NInf, KInf*NInf); 
pInf = NInf/(p+1)-1;  qOrdInf = max(2*pInf+1, 1); [~,WInf] = quadRuleInf1D(qOrdInf, beta);
for nn = 1 : 3
  for np = 1 : 3
    for r = 1 : length(WInf)
      if nn ~= 2 && np ~=2
        HLPow1o2 = gCdiagInf{1,nn}(:,r).^0.5;
        HRPow1o2 = gCoffdiagInf{1,nn,np}(:,r).^0.5;
        lambda =   abs( ( ( HRPow1o2 .* gCdiagInf{2,nn}(:,r) + HLPow1o2 .* gCoffdiagInf{2,nn,np}(:,r) ) .* gInf.nuE0Inf(:,nn,1)   ...
                        + ( HRPow1o2 .* gCdiagInf{3,nn}(:,r) + HLPow1o2 .* gCoffdiagInf{3,nn,np}(:,r) ) .* gInf.nuE0Inf(:,nn,2) ) ...
                        ./( HRPow1o2 .* gCdiagInf{1,nn}(:,r) + HLPow1o2 .* gCoffdiagInf{1,nn,np}(:,r)) );
        lambda = ( min(lambda, 0) + max(lambda, 0) ...
                 + sqrt( gConst * (gCdiagInf{1,nn}(:,r).^1.5 + gCoffdiagInf{1,nn,np}(:,r).^1.5) ./ (HLPow1o2 + HRPow1o2) ) ) ...
                 .* gInf.detE0I(:,nn);
        globV = globV + 0.5 * kron( spdiags( gInf.markE0Iaux{nn,np} .* lambda, 0, KInf, KInf ), PhiPhidiag(:,:,r,nn) ) ...
                      - 0.5 * kron( bsxfun(@times, gInf.markE0IE0I{nn,np}, gInf.markE0Iaux{nn,np} .* lambda ), PhiPhioffdiag(:,:,r,nn,np) );
      end % if
    end % for
  end % for
end % for
end % function
