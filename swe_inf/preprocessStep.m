
function pd = preprocessStep(pd, nStep)
pd = computeSolutionOnQuad(pd);
pd = computeSolutionOnQuadInf(pd);
cDGInf = pd.cDiscInf;
KInf = pd.KInf;
NInf = pd.NInf;
gInf = pd.gInf;
beta = pd.beta;
p = pd.p;
pInf = pd.pInf;
gConst = pd.gConst;
tauConst = pd.bottomFrictionCoef;
t = pd.t;
%% coefficients of primary unknowns at previous time step
pd.sysY =      [reshape(cDGInf(:,:,1).', KInf*NInf, 1); reshape(cDGInf(:,:,2).', KInf*NInf, 1); reshape(cDGInf(:,:,3).', KInf*NInf, 1)  ];
%% assembly of time-dependent global area integral ontributions of non-linearity, discretized explicitly
pd.globEInf    = assembleGlobEInf   (gInf, pd.PhiPhi2DInf,                  cDGInf, tauConst, p, beta, pd.solutionOnQuadInf.gC2DInf);
pd.globFInf    = assembleGlobFInf   (gInf, pd.hatGInf,                      cDGInf,   gConst, p, beta);
pd.globFExpInf = assembleGlobFExpInf(gInf, pd.hatDInf,                      cDGInf,   gConst, p, beta);
pd.globOInf    = assembleGlobOInf   (gInf, pd.gradPhiPhi2DInf,              cDGInf,           p, beta, pd.solutionOnQuadInf.gC2DInf);
pd.globOExpInf = assembleGlobOExpInf(gInf, pd.PhiPhi2DInf,                  cDGInf,           p, beta, pd.solutionOnQuadInf.gC2DInf);
%% assembly of time-dependent global interior edge integral contributions of non-linearity, discretized explicitly
pd.globRInf = assembleGlobRInf(gInf, pd.hatRdiagInf,   pd.hatRoffdiagInf,   cDGInf,   gConst         );
pd.globUInf = assembleGlobUInf(gInf, pd.PhiPhidiagInf, pd.PhiPhioffdiagInf, cDGInf          , p, beta, pd.solutionOnQuadInf.gCdiagInf, pd.solutionOnQuadInf.gCoffdiagInf);
pd.globVInf = assembleGlobVInf(gInf, pd.PhiPhidiagInf, pd.PhiPhioffdiagInf, cDGInf,   gConst, p, beta, pd.solutionOnQuadInf.gCdiagInf, pd.solutionOnQuadInf.gCoffdiagInf);
%% assembly of time-dependent global land boundary edge integral contributions of non-linearity, discretized explicitly
pd.globRLInf  = assembleGlobRLInf(gInf, pd.hatRdiagInf, cDGInf, gConst);
%% assembly of time-dependent global open sea boundary edge integral contributions of non-linearity, discretized explicitly
pd.globUOSInf = assembleGlobUOSInf(gInf, @(x1,x2) pd.uCont(x1,x2,t), @(x1,x2) pd.vCont(x1,x2,t), @(x1,x2) pd.HOSAlg(x1,x2,t), NInf,p, beta, pd.basesOnQuadInf);%neuer Zeitschritt?!
pd.globROSInf = assembleGlobROSInf(gInf, @(x1,x2) pd.HOSAlg(x1,x2,t), NInf, p, gConst, beta, pd.basesOnQuadInf);
pd.globLInf    = cell(3,1);
F1DGInf     = projectAlg2Inf(gInf, @(x1,x2) pd.f1Cont(x1,x2,t), 2*p, 2*pInf, pd.hatMInf, beta, pd.basesOnQuadInf);
F2DGInf     = projectAlg2Inf(gInf, @(x1,x2) pd.f2Cont(x1,x2,t), 2*p, 2*pInf, pd.hatMInf, beta, pd.basesOnQuadInf);
pd.globLInf{1} = pd.globMInf * reshape(F1DGInf', KInf*NInf, 1);
pd.globLInf{2} = pd.globMInf * reshape(F2DGInf', KInf*NInf, 1);
pd.globLInf{3} = zeros(KInf*NInf,1);
if pd.isRhsAvail
  F0DGInf = projectAlg2Inf(gInf, @(x1,x2) pd.f0Cont(x1,x2,t), 2*p, 2*pInf, pd.hatMInf, beta, pd.basesOnQuadInf);
  pd.globLInf{3} = pd.globMInf * reshape(F0DGInf', KInf*NInf, 1);
end % if
[pd.globROSRiemInf, pd.globUOSRiemInf, pd.globVOSRiemInf] = additionalOSRiemInf(gInf, pd.HOSAlg, pd.uCont, pd.vCont, pd.g.markE0TOS, pd.hatRdiagInf, pd.PhiPhidiagInf, cDGInf, p, gConst, t, beta, pd.solutionOnQuadInf.gCdiagInf, pd.solutionOnQuadInf.gCoffdiagInf, pd.basesOnQuadInf);
end % function
