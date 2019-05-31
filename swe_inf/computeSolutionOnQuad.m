function pd = computeSolutionOnQuad(pd)
%% computes values of unknowns on all necessary quadrature points at previous time step
g =pd.g;
cDG = pd.cDisc;
K = g.numT; N = size(cDG, 2);
p = (sqrt(8*N+1)-3)/2;
if p > 0, qOrd = 2*p; else qOrd = 1; end
[~, ~, W] = quadRule2D(qOrd);
R = length(W);
pd.solutionOnQuad.gC2D = cell(3,1);
pd.solutionOnQuad.gC2D{1} = zeros(K,R);
pd.solutionOnQuad.gC2D{2} = zeros(K,R);
pd.solutionOnQuad.gC2D{3} = zeros(K,R);
for l = 1:N
  pd.solutionOnQuad.gC2D{1} = pd.solutionOnQuad.gC2D{1} + cDG(:,l,1) * pd.basesOnQuad.phi2D{qOrd}(:,l).';
  pd.solutionOnQuad.gC2D{2} = pd.solutionOnQuad.gC2D{2} + cDG(:,l,2) * pd.basesOnQuad.phi2D{qOrd}(:,l).';
  pd.solutionOnQuad.gC2D{3} = pd.solutionOnQuad.gC2D{3} + cDG(:,l,3) * pd.basesOnQuad.phi2D{qOrd}(:,l).';
end % for
qOrd = 2*p+1;
[~, W] = quadRule1D(qOrd);
R = length(W);
pd.solutionOnQuad.gCdiag = cell(3,3);
for n = 1 : 3
  pd.solutionOnQuad.gCdiag{1,n} = zeros(K,R);
  pd.solutionOnQuad.gCdiag{2,n} = zeros(K,R);
  pd.solutionOnQuad.gCdiag{3,n} = zeros(K,R);
  for l = 1 : N
    pd.solutionOnQuad.gCdiag{1,n} = pd.solutionOnQuad.gCdiag{1,n} + cDG(:,l,1) * pd.basesOnQuad.phi1D{qOrd}(:,l,n).';
    pd.solutionOnQuad.gCdiag{2,n} = pd.solutionOnQuad.gCdiag{2,n} + cDG(:,l,2) * pd.basesOnQuad.phi1D{qOrd}(:,l,n).';
    pd.solutionOnQuad.gCdiag{3,n} = pd.solutionOnQuad.gCdiag{3,n} + cDG(:,l,3) * pd.basesOnQuad.phi1D{qOrd}(:,l,n).';
  end % for
end % for
pd.solutionOnQuad.gCoffdiag = cell(3,3,3);
for nn = 1 : 3
  for np = 1 : 3
    pd.solutionOnQuad.gCoffdiag{1,nn,np} = zeros(K,R);
    pd.solutionOnQuad.gCoffdiag{2,nn,np} = zeros(K,R);
    pd.solutionOnQuad.gCoffdiag{3,nn,np} = zeros(K,R);
    for l = 1 : N
      pd.solutionOnQuad.gCoffdiag{1,nn,np} = pd.solutionOnQuad.gCoffdiag{1,nn,np} + g.markE0TE0T{nn,np} * cDG(:,l,1) * pd.basesOnQuad.thetaPhi1D{qOrd}(:,l,nn,np).';
      pd.solutionOnQuad.gCoffdiag{2,nn,np} = pd.solutionOnQuad.gCoffdiag{2,nn,np} + g.markE0TE0T{nn,np} * cDG(:,l,2) * pd.basesOnQuad.thetaPhi1D{qOrd}(:,l,nn,np).';
      pd.solutionOnQuad.gCoffdiag{3,nn,np} = pd.solutionOnQuad.gCoffdiag{3,nn,np} + g.markE0TE0T{nn,np} * cDG(:,l,3) * pd.basesOnQuad.thetaPhi1D{qOrd}(:,l,nn,np).';
    end % for
  end % for
end % for
end % function
