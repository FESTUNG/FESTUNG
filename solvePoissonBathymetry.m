function solvePoissonBathymetry(g, diffusionType, depth)

N = 3;

bdrInd = [1 4 78:102];
% bdrInd = 1:102;

rhs = zeros(g.numV,1);
rhs(bdrInd) = -depth(bdrInd,end);

if ~strcmp(diffusionType, 'TVD')
  [A, step] = getDiffusionMatrix(g, diffusionType, N);
  for i = 1: length(bdrInd)
    A(bdrInd(i),bdrInd(i)) = 1;
    A(bdrInd(i),1:bdrInd(i)-1) = 0;
    A(bdrInd(i),bdrInd(i)+1:end) = 0;
  end
  
  b=A\rhs;
else
  b = zeros(g.numV,1);
  changeL2 = Inf; iter = 0;
%   refElemDphiLagrDphiLagr = integrateRefElemDphiLagrDphiLagr(N);
%   AP = assembleMatElemDphiLagrDphiLagr(g, refElemDphiLagrDphiLagr);
  
  while changeL2 > 1e-10 && iter < Inf
    [A, step] = getDiffusionMatrix(g, diffusionType, N, b);
%     diffKoeff = @(x) (x(1)^2+x(2)^2+1e-5)^-0.5;
%     refElemDphiLagrDphiLagr = integrateRefElemDphiLagrDphiLagr(N);
%     ATVD = assembleMatElemDphiLagrDphiLagrVal(g, refElemDphiLagrDphiLagr, b, diffKoeff);
    for i = 1: length(bdrInd)
      A(bdrInd(i),bdrInd(i)) = 1;
      A(bdrInd(i),1:bdrInd(i)-1) = 0;
      A(bdrInd(i),bdrInd(i)+1:end) = 0;
      
%       ATVD(bdrInd(i),bdrInd(i)) = 1;
%       ATVD(bdrInd(i),1:bdrInd(i)-1) = 0;
%       ATVD(bdrInd(i),bdrInd(i)+1:end) = 0;
    end
    bOld = b;
%     b = b + ATVD \ (rhs - A * b);
    b = A \ rhs;
    changeL2 = norm(b-bOld);
    iter = iter + 1;
  end
  if changeL2 > 1e-10
    fprintf('Fixed point method does not converge, error %1.2e.\n', changeL2)
  else
    iter
  end
end

[Q1, Q2] = quadRule2D(2); 
basesOnQuad = computeBasesOnQuad(N, struct);
refElemPhiPhi = integrateRefElemPhiPhi(N, basesOnQuad);
bDisc = projectDataQ0T2DataDisc(b(g.V0T) * [1-Q1-Q2; Q1; Q2], 2, refElemPhiPhi, basesOnQuad);
bLagr = projectDataDisc2DataLagr(bDisc);

if length(bdrInd) == 102
  step = step+3;
end
visualizeDataLagr(g,bLagr,'lakeAtRest','lakeAtRest',step,'vtk')

end

function [A, step] = getDiffusionMatrix(g, diffusionType, N, b)

switch diffusionType
  case 'artificial'
    refElemPhiLagrPhiLagr = integrateRefElemPhiLagrPhiLagr();
    massMatLumped = assembleVecElemPhiLagr(g);
    A = -(assembleMatElemPhiLagrPhiLagr(g, refElemPhiLagrPhiLagr) - spdiags(massMatLumped, 0, g.numV, g.numV));
    step = 0;
  case 'poisson'
    refElemDphiLagrDphiLagr = integrateRefElemDphiLagrDphiLagr(N);
    A = assembleMatElemDphiLagrDphiLagr(g, refElemDphiLagrDphiLagr);
    step = 1;
  case 'TVD'
    diffKoeff = @(x) (x(1)^2+x(2)^2+1e-10)^-0.5;
    refElemDphiLagrDphiLagr = integrateRefElemDphiLagrDphiLagr(N);
    A = assembleMatElemDphiLagrDphiLagrVal(g, refElemDphiLagrDphiLagr, b, diffKoeff);
    step = 2;
end
end