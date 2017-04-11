function problemData = assembleJumpTerms(problemData, nSubStep)
t = problemData.t(nSubStep);
K = problemData.g.numT;
barK = problemData.g.g1D.numT;

[Q,~] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);
hDCont = @(x1) problemData.hDCont(t,x1);
u1DCont = @(x1,x2) problemData.u1DCont(t,x1,x2);

u1Q0E0Tint = { zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % cDisc{2} in quad points of edges
u1Q0E0TE0T = { zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % cDisc{2} of neighboring element in quad points of edges
u1Q0E0TbdrRiem = { zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % u1D in quad points of boundary edges with Riemann solver
for n = 1 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,n) * problemData.cDiscRK{nSubStep, 2}.', K * numQuad1D, 1);
  [Q1, Q2] = gammaMapTetra(n, Q);
  X1 = problemData.g.mapRef2Phy(1, Q1, Q2);
  X2 = problemData.g.mapRef2Phy(2, Q1, Q2);
  u1Q0E0TbdrRiem{n}(problemData.g.markE0TbdrRiem(:, n)) = u1DCont(X1(problemData.g.markE0TbdrRiem(:, n)), X2(problemData.g.markE0TbdrRiem(:, n)));
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,mapLocalEdgeTetra(n)) * problemData.cDiscRK{nSubStep, 2}.';
  u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', K * numQuad1D, 1);
end % for nn

problemData.globKu = zeros(K * problemData.N, 1);
problemData.globKh = zeros(K * problemData.N, 1);
problemData.barGlobKh = zeros(barK * problemData.barN, 1);
for n = 3 : 4
  nn1D = 5 - n; np1D = 5 - mapLocalEdgeTetra(n);
  
  [i, j] = find(problemData.g.g1D.markV0TbdrRiem(:, nn1D));
  hV0T1DbdrRiem = sparse(i, j, hDCont(problemData.g.g1D.coordV0T(problemData.g.g1D.markV0TbdrRiem(:, nn1D), nn1D, 1)), barK, 1);
  
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,nn1D) + problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) + hV0T1DbdrRiem);
  hJmpE0T = problemData.g.g1D.markT2DT * ( ( problemData.hV0T1D(:,nn1D) - problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) - hV0T1DbdrRiem ) ./ problemData.hSmoothV0T1D(:,nn1D) );
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n} + u1Q0E0TbdrRiem{n});
  lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  hJmpLambdaE0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D, 1));
    
  problemData.globKu = problemData.globKu + problemData.globSu{n} * ( lambdaQ0E0T .* (u1Q0E0Tint{n} - u1Q0E0TE0T{n} - u1Q0E0TbdrRiem{n}) );
  problemData.globKh = problemData.globKh + problemData.globSh{n} * hJmpLambdaE0T;
  problemData.barGlobKh = problemData.barGlobKh + problemData.barGlobS{n} * hJmpLambdaE0T;
end % for n
end % function