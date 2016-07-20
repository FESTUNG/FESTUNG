function pd = initializeProblem(pd)
p = pd.p;
K = pd.K;
N = pd.N;
pInf = pd.pInf;
KInf = pd.KInf;
NInf = pd.NInf;

%% Primary unknowns H, uH, vH (each K*N).
if pd.isSolutionAvail
  xi0Cont = @(x1,x2) pd.xiCont(x1,x2,pd.t0);
  u0Cont = @(x1,x2) pd.uCont(x1,x2,pd.t0);
  v0Cont = @(x1,x2) pd.vCont(x1,x2,pd.t0);
end % if

pd.cDiscInf = initializeCDGInf(pd.gInf,@(x1,x2) xi0Cont(x1,x2),@(x1,x2) u0Cont(x1,x2),@(x1,x2) v0Cont(x1,x2), pd.zbInfQuad, 2*p, 2*pInf, pd.hatMInf, pd.beta, pd.basesOnQuadInf);
pd.cDisc = zeros(K,N,3); %braucht man eigentlich nicht mehr, habe es aber mal dringelassen, damit ihr wisst wo es vorkommt (ich lass auch bei den ehemals globalen Variablen mal die offdiag-Strukturen, ihr koennt sie ja rausnehmen)
pd.cDisc(:,:,1) = projectAlg2DG(pd.g, @(x1,x2)  xi0Cont(x1,x2) - pd.zbCont(x1,x2)                  , 2*p, pd.hatM, pd.basesOnQuad);
pd.cDisc(:,:,2) = projectAlg2DG(pd.g, @(x1,x2) (xi0Cont(x1,x2) - pd.zbCont(x1,x2)) .* u0Cont(x1,x2), 2*p, pd.hatM, pd.basesOnQuad);
pd.cDisc(:,:,3) = projectAlg2DG(pd.g, @(x1,x2) (xi0Cont(x1,x2) - pd.zbCont(x1,x2)) .* v0Cont(x1,x2), 2*p, pd.hatM, pd.basesOnQuad);

pd.cDisc =    correctHeight(p,  pd.cDisc, pd.minTol, pd.sysMinValueCorrection, 0);
pd.cDiscInf = correctHeightInf(p, pInf, pd.cDiscInf, pd.minTol, pd.sysMinValueCorrectionInf, 0, pd.beta);

%% Initialize waitbar.
if pd.isWaitbar
  str  = strcat( ['% done. Simulating refinement level =', ' ', num2str(pd.refinement), ', p =', ' ', num2str(p), ', pInf =', ' ', num2str(pInf), '.' ]);
  pd.waitbar = waitbar(0, strcat([ 'Time stepping:', ' ', num2str(0), str ]));
end % if

%% Initialize time stepping
pd.isFinished = false;
pd.t = pd.t0;
end % function