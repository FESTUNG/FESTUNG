%%%% script for testing advection equation
clc; clear
g = domainSquare(2^-5);
N = 3; p = (sqrt(8*N+1)-3)/2
qOrd = max(2*p,1);
b = computeBasesOnQuad(N, struct);
u = projectFuncCont2DataDisc(g,@(x,y) 0.5-y,qOrd,eye(N),b);
v = projectFuncCont2DataDisc(g,@(x,y) x-0.5,qOrd,eye(N),b);

G = @(x1, x2, x1_0, x2_0) (1/0.15) * sqrt((x1-x1_0).^2 + (x2-x2_0).^2);
z = @(x1, x2) ((x1 - 0.5).^2 + (x2 - 0.75).^2 <= 0.0225 & (x1 <= 0.475 | x1 >= 0.525 | x2 >= 0.85)) + ...
                    (1-G(x1, x2, 0.5, 0.25)) .* ((x1 - 0.5).^2 + (x2 - 0.25).^2 <= 0.0225) + ...
                    0.25*(1+cos(pi*G(x1, x2, 0.25, 0.5))).*((x1 - 0.25).^2 + (x2 - 0.5).^2 <= 0.0225);

zLagr = zeros(g.numV,1);
for i = 1:g.numV
  zLagr = z(g.coordV(:,1),g.coordV(:,2));
end

MLumpedInv = assembleVecElemPhiLagr(g).^-1;
refElemDphiLagrPhi = integrateRefElemDphiLagrPhi(N, b);
A = assembleMatElemDphiLagrPhi(g,refElemDphiLagrPhi);

[Q1, Q2] = quadRule2D(qOrd);

t= 0; dt = 0.000002; outp = 0; outFreq = 100000; outNum = 0;
while t < 2*pi
  
  zV0T = zLagr(g.V0T);
  dataQ0T = (u * b.phi2D{qOrd}') .* (zV0T * [1 - Q1 - Q2; Q1; Q2]);
  U = reshape(projectDataQ0T2DataDisc(dataQ0T, qOrd, eye(N),b), [], 1);
  dataQ0T = (v * b.phi2D{qOrd}') .* (zV0T * [1 - Q1 - Q2; Q1; Q2]);
  V = reshape(projectDataQ0T2DataDisc(dataQ0T, qOrd, eye(N),b), [], 1);
  
%   if p == 0
%     U = sum(zV0T,2)/3;
%   else
%     U = reshape(zV0T, [], 1); 
%   end
%   
%   V = U;
%   U = zeros(size(V));
  
  if mod(outp, outFreq) == 0
    visualizeDataLagr(g, zV0T,'rotationCG','rotationCG', outNum,'vtk');
    outNum = outNum + 1;
  end
  
  zLagr = zLagr + dt * MLumpedInv .* (A{1} * U + A{2} * V);
  
  t = t + dt; outp = outp+1;
end