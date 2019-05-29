
function pd = preprocessProblem(pd)
%% Triangulation
switch pd.gridSource

    
  case 'hierarchical'
    X1 = [0 100 100 0]; X2 = [0 0 100 100];
    pd.g = domainHierarchy(X1, X2, pd.hmax, pd.refinement);
    
    % Set edge types
%     pd.g.idE = zeros(pd.g.numE,1);
%     pd.g.idE(pd.g.baryE(:, 2) == 0) = 1; % south
%     pd.g.idE(pd.g.baryE(:, 1) == 1) = 4; % east
%     pd.g.idE(pd.g.baryE(:, 2) == 1) = 1; % north
%     pd.g.idE(pd.g.baryE(:, 1) == 0) = 1; % west
%     pd.g.idE0T = pd.g.idE(pd.g.E0T);
    pd.g.markE0TL = pd.g.idE0T == 1 | pd.g.idE0T == 3 | pd.g.idE0T == 4;
    
    % Store edge counts
    pd.g.numEint = sum(pd.g.idE == 0);
    pd.g.numEbdrL = sum(pd.g.idE == 1);
    pd.g.numEbdrRA = sum(pd.g.idE == 2);
    pd.g.numEbdrRI = sum(pd.g.idE == 3);
    pd.g.numEbdrOS = sum(pd.g.idE == 4);
    
    pd.g.markE0Taux = cell(3,3); % auxiliary cell of vectors of length K, needed in some routines
    for nn = 1:3
      for np = 1:3
        pd.g.markE0Taux{nn,np} = pd.g.markE0TE0T{nn,np} * ones(pd.g.numT,1);
      end
    end
    
    pd.g.markE0Tint = pd.g.idE0T == 0;                   % [K x 3] mark local edges that are interior
    pd.g.markE0TOS  = ~(pd.g.markE0Tint | pd.g.markE0TL);% [K x 3] mark local edges on the open sea boundary
    
    pd.gInf  = gridInf(pd);
    
  otherwise
    error('Invalid gridSource given.')
end % switch

%% Globally constant parameters
pd.K = pd.g.numT; % number of triangles
pd.N = nchoosek(pd.p + 2, pd.p); % number of local DOFs
K = pd.K;
N = pd.N;

pd.KInf = pd.gInf.numInfElem;
pd.NInf = (pd.p+1)*(pd.pInf+1);
KInf = pd.KInf;
NInf = pd.NInf;

%% lookup table for basis function
pd.basesOnQuadInf = computeBasesOnQuadInf(N, NInf, pd);
%% lookup table for basis function
pd.basesOnQuad = computeBasesOnQuad(N, struct, [max(2*pd.p, 1), 2*pd.p+1]);
%% compute gInf and zbLin
gInf  = gridInf(pd);
KInf  = gInf.numInfElem;
zbLin = projectZb2Lin(pd.gInf, pd.zbOSAlg);%<---Koeffizienten fuer phiLin auf Ref-Element
%% computation of matrices on the reference inf-element
pd.hatDInf          = computeHatDInf         (NInf, pd);
pd.hatGInf          = computeHatGInf         (NInf, pd);
hatGzbInf        = computeHatGzbInf       (NInf, pd);
hatHInf          = computeHatHInf         (NInf, pd);
pd.hatMInf          = computeHatMInf         (NInf, pd);
pd.hatRdiagInf      = computeHatRdiagInf     (NInf, pd);
pd.hatRoffdiagInf   = computeHatRoffdiagInf  (NInf, pd);
hatSdiagInf      = computeHatSdiagInf     (NInf, pd);
hatSoffdiagInf   = computeHatSoffdiagInf  (NInf, pd);
pd.PhiPhi2DInf      = computePhiPhiInf2D     (NInf, pd);
pd.gradPhiPhi2DInf  = computeGradPhiPhiInf2D (NInf, pd);
pd.PhiPhidiagInf    = computePhiPhidiagInf   (NInf, pd);
pd.PhiPhioffdiagInf = computePhiPhioffdiagInf(NInf, pd);
%% computation of matrices on the reference triangle
pd.hatM = integrateRefElemPhiPhi(N, pd.basesOnQuad);
%% L2 projections of time-independent algebraic coefficients
fcInf = projectAlg2Inf(gInf, @(x1,x2) pd.fcCont(x1,x2), 2*pd.p, 2*pd.pInf, pd.hatMInf, pd.beta, pd.basesOnQuadInf);
% only use zbDG/zbInfQuad to compute xi from c1 and for error analysis
pd.zbInfQuad = zbAtInfQuadPoints(pd.gInf, pd.zbOSAlg, pd.p, pd.pInf, pd.beta);%constant extension!



%% System matrix for correction of min value exceedence.
pd.sysMinValueCorrection = [ phi(1,0,0) phi(1,1,0) phi(1,0,1) ; ...
                             phi(2,0,0) phi(2,1,0) phi(2,0,1) ; ...
                             phi(3,0,0) phi(3,1,0) phi(3,0,1) ];
                         
pd.sysMinValueCorrectionInf = [ phiInf(1,1,0,0,pd.beta) phiInf(1,1,0,1,pd.beta) ; ...
                                phiInf(2,1,0,0,pd.beta) phiInf(2,1,0,1,pd.beta) ];



%% Assembly of time-independent global matrices
% Element matrices
pd.globDInf    = assembleGlobDInf    (pd.gInf,             pd.hatDInf,              fcInf           , pd.p, pd.beta);
pd.globGInf    = assembleGlobGInf    (pd.gInf,             hatGzbInf,               zbLin, pd.gConst, pd.p, pd.beta);
pd.globHInf    = assembleGlobHInf    (pd.gInf,             hatHInf,                                   pd.p, pd.beta);
pd.globHExpInf = assembleGlobHExpInf (pd.gInf,             pd.hatMInf,                                pd.p, pd.beta);
pd.globMInf    = assembleGlobMInf    (pd.gInf,             pd.hatMInf                               , pd.p, pd.beta);
% Boundary edge matrices
pd.globQInf    = assembleGlobQInf    (pd.gInf, hatSdiagInf, hatSoffdiagInf);
pd.globQOSInf  = assembleGlobQOSInf  (pd.gInf, hatSdiagInf);

% Derived system matrices
pd.sysW = blkdiag(pd.globMInf, pd.globMInf, pd.globMInf);

end % function
