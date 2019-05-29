function pd = computeSolutionOnQuadInf(pd)

gInf = pd.gInf; 
cDGInf = pd.cDiscInf;
cDG = pd.cDisc;
K   = pd.K;
N = pd.N;
beta = pd.beta;
%% computes values of unknowns on all necessary quadrature points at previous time step
KInf = gInf.numInfElem; NInf = size(cDGInf, 2); 
p = (sqrt(8*N+1)-3)/2;       
pInf = NInf/(p+1)-1;
if p > 0,    qOrd    = 2*p;    else qOrd    = 1; end
if pInf > 0, qOrdInf = 2*pInf; else qOrdInf = 1; end
[~, ~, W] = quadRuleInf2D(qOrd, qOrdInf, beta);
R = length(W);
pd.solutionOnQuadInf.gC2DInf = cell(3,1);
pd.solutionOnQuadInf.gC2DInf{1} = zeros(KInf,R);
pd.solutionOnQuadInf.gC2DInf{2} = zeros(KInf,R);
pd.solutionOnQuadInf.gC2DInf{3} = zeros(KInf,R);
for l = 1:NInf
  pd.solutionOnQuadInf.gC2DInf{1}(:,:) = pd.solutionOnQuadInf.gC2DInf{1}(:,:) + cDGInf(:,l,1) * pd.basesOnQuadInf.gPhiInf2D{qOrdInf}(:,l).';
  pd.solutionOnQuadInf.gC2DInf{2}(:,:) = pd.solutionOnQuadInf.gC2DInf{2}(:,:) + cDGInf(:,l,2) * pd.basesOnQuadInf.gPhiInf2D{qOrdInf}(:,l).';
  pd.solutionOnQuadInf.gC2DInf{3}(:,:) = pd.solutionOnQuadInf.gC2DInf{3}(:,:) + cDGInf(:,l,3) * pd.basesOnQuadInf.gPhiInf2D{qOrdInf}(:,l).';
end % for

% qOrd    = 2*p+1;
qOrdInf = 2*pInf+1;
[~, W] = quadRuleInf1D(qOrdInf, beta);
R = length(W);
pd.solutionOnQuadInf.gCdiagInf = cell(3,3);
for n = 1 : 3
%   if n == 2
%     [~, W] = quadRule1D(qOrd);
%   else
%     [~, W] = quadRuleInf1D(qOrdInf, beta);
%   end % if
%   R = length(W);
  pd.solutionOnQuadInf.gCdiagInf{1,n} = zeros(KInf,R);
  pd.solutionOnQuadInf.gCdiagInf{2,n} = zeros(KInf,R);
  pd.solutionOnQuadInf.gCdiagInf{3,n} = zeros(KInf,R);
  for l = 1 : NInf
    pd.solutionOnQuadInf.gCdiagInf{1,n}(:,:) = pd.solutionOnQuadInf.gCdiagInf{1,n}(:,:) + cDGInf(:,l,1) * pd.basesOnQuadInf.gPhiInf1D{qOrdInf}(:,l,n).';
    pd.solutionOnQuadInf.gCdiagInf{2,n}(:,:) = pd.solutionOnQuadInf.gCdiagInf{2,n}(:,:) + cDGInf(:,l,2) * pd.basesOnQuadInf.gPhiInf1D{qOrdInf}(:,l,n).';
    pd.solutionOnQuadInf.gCdiagInf{3,n}(:,:) = pd.solutionOnQuadInf.gCdiagInf{3,n}(:,:) + cDGInf(:,l,3) * pd.basesOnQuadInf.gPhiInf1D{qOrdInf}(:,l,n).';
  end % for
end % for
pd.solutionOnQuadInf.gCoffdiagInf = cell(3,3,3);
for nn = 1 : 3
  for np = 1 : 3
    if nn == 2
      pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np} = zeros(KInf,R);
      pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np} = zeros(KInf,R);
      pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np} = zeros(KInf,R);
      for l = 1 : N
        pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np}(:,:) + gInf.markE0TE0I{nn,np} * cDG(:,l,1) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
        pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np}(:,:) + gInf.markE0TE0I{nn,np} * cDG(:,l,2) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
        pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np}(:,:) + gInf.markE0TE0I{nn,np} * cDG(:,l,3) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
      end % for
    elseif np == 2
      pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np} = zeros(K,R);
      pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np} = zeros(K,R);
      pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np} = zeros(K,R);
      for l = 1 : NInf
        pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np}(:,:) + gInf.markE0TE0I{nn,np} * cDGInf(:,l,1) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
        pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np}(:,:) + gInf.markE0TE0I{nn,np} * cDGInf(:,l,2) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
        pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np}(:,:) + gInf.markE0TE0I{nn,np} * cDGInf(:,l,3) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
      end % for
    else
      pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np} = zeros(KInf,R);
      pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np} = zeros(KInf,R);
      pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np} = zeros(KInf,R);
      for l = 1 : NInf
        pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{1,nn,np}(:,:) + gInf.markE0IE0I{nn,np} * cDGInf(:,l,1) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
        pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{2,nn,np}(:,:) + gInf.markE0IE0I{nn,np} * cDGInf(:,l,2) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
        pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np}(:,:) = pd.solutionOnQuadInf.gCoffdiagInf{3,nn,np}(:,:) + gInf.markE0IE0I{nn,np} * cDGInf(:,l,3) * pd.basesOnQuadInf.gThetaPhiInf1D{qOrdInf}(:,l,nn,np).';
      end % for
    end % if
  end % for
end % for

pd.solutionOnQuadInf.gCoffdiagT2I2 = cell(3);%nn=np=2, nn index of triangle, np index of inf
pd.solutionOnQuadInf.gCoffdiagT2I2{1} = zeros(K,R);
pd.solutionOnQuadInf.gCoffdiagT2I2{2} = zeros(K,R);
pd.solutionOnQuadInf.gCoffdiagT2I2{3} = zeros(K,R);
for l = 1 : NInf
  pd.solutionOnQuadInf.gCoffdiagT2I2{1}(:,:) = pd.solutionOnQuadInf.gCoffdiagT2I2{1}(:,:) + gInf.markE0TE0I{2,2}' * cDGInf(:,l,1) * pd.basesOnQuadInf.gThetaPhiLinkT2I21D{qOrdInf}(:,l).';
  pd.solutionOnQuadInf.gCoffdiagT2I2{2}(:,:) = pd.solutionOnQuadInf.gCoffdiagT2I2{2}(:,:) + gInf.markE0TE0I{2,2}' * cDGInf(:,l,2) * pd.basesOnQuadInf.gThetaPhiLinkT2I21D{qOrdInf}(:,l).';
  pd.solutionOnQuadInf.gCoffdiagT2I2{3}(:,:) = pd.solutionOnQuadInf.gCoffdiagT2I2{3}(:,:) + gInf.markE0TE0I{2,2}' * cDGInf(:,l,3) * pd.basesOnQuadInf.gThetaPhiLinkT2I21D{qOrdInf}(:,l).';
end % for
end % function
