function problemData = preprocessStep(problemData, nStep)
t = nStep * problemData.tau;
barK = problemData.g.g1D.numT;

%% Apply mesh adaptation to free surface movement
problemData = adaptMesh(problemData);

%% L2-projections of algebraic coefficients and right hand side.
DDisc = cellfun(@(c) execin('darcyVert/projectFuncCont2DataDiscTrap', problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
        problemData.hatM{1}, problemData.basesOnQuad2D), problemData.DCont, 'UniformOutput', false);
problemData.globLu = problemData.globM * ...
        reshape(execin('darcyVert/projectFuncCont2DataDiscTrap', problemData.g, @(x1,x2) problemData.fuCont(t,x1,x2), ...
        problemData.qOrd, problemData.hatM{1}, problemData.basesOnQuad2D).', [], 1);
problemData.globLh = problemData.barGlobM * reshape(projectFuncCont2DataDisc1D(problemData.g.g1D, @(x1) problemData.fhCont(t,x1), ...
        problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D).', [], 1);
%% Determine quadrature rules
[Q,W] = quadRule1D(problemData.qOrd); 
numQuad1D = length(W);
mapE0E = [2 1 4 3];
      
%% Create lookup tables for solution on quadrature points
heightV0T1D = [ problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), 4, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 2), ...
                problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), 3, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 2) ];
heightQ0T1D = heightV0T1D(:,1) * (1-Q) + heightV0T1D(:,2) * Q;
              
u1Q0E0Tint = cell(4,1); % cDisc{2} in quad points of edges
u1Q0E0TE0T = cell(4,1); % cDisc{2} of neighboring element in quad points of edges
for n = 1 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D(:,:,n) * problemData.cDisc{2}.', problemData.g.numT * numQuad1D, 1);
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D(:,:,mapE0E(n)) * problemData.cDisc{2}.';
  u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', problemData.g.numT * numQuad1D, 1);
end % for nn

%% Compute depth averaged velocity
barU1Disc = { zeros(problemData.g.g1D.numT, problemData.barN), zeros(problemData.g.g1D.numT, problemData.barN) };
for s = 1 : 2
  for j = 1 : problemData.barN
    i = execin('darcyVert/mapTensorProductIndex', j, 1);
    barU1Disc{s}(:, j) = problemData.g.g1D.markT2DT.' * (problemData.cDisc{2}(:,i) .* problemData.g.J0T{s}(:,2,2));
  end % for j
end % for s

%% Assembly of time-dependent global matrices.
problemData.globE = execin('darcyVert/assembleMatElemTrapDphiPhiFuncDisc', problemData.g, problemData.hatG, problemData.cDisc(2:3));
problemData.globE = problemData.globE{1} + problemData.globE{2};
problemData.globG = execin('darcyVert/assembleMatElemTrapDphiPhiFuncDisc', problemData.g, problemData.hatG, DDisc);

problemData.globR = execin('darcyVert/assembleMatEdgeTrapPhiPhiFuncDiscNu', problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);
problemData.globP = execin('darcyVert/assembleMatEdgeTrapPhiPhiFuncDiscNu', problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDisc(2:3));
problemData.globP = problemData.globP{1} + problemData.globP{2};

problemData.globJu = zeros(problemData.g.numT * problemData.N, 1);
problemData.globJh = zeros(problemData.g.numT * problemData.N, 1);
problemData.barGlobJh = zeros(problemData.g.g1D.numT * problemData.barN, 1);
for n = 3 : 4
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,n-2) + spdiags(ones(barK, 1), 7-2*n, barK, barK) * problemData.hV0T1D(:,mapE0E(n)-2) );
  hJmpE0T = problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,n-2) - spdiags(ones(barK, 1), 7-2*n, barK, barK) * problemData.hV0T1D(:,mapE0E(n)-2) );
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n});
  lambdaE0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  lambdaE0T = 0 * lambdaE0T;
  hJmpLambdaE0T = lambdaE0T .* kron(hJmpE0T, ones(numQuad1D,1));
    
  problemData.globJu = problemData.globJu + problemData.globS{n} * ( lambdaE0T .* (u1Q0E0Tint{n} - u1Q0E0TE0T{n}) );
  problemData.globJh = problemData.globJh + problemData.globS{n} * hJmpLambdaE0T;
  problemData.barGlobJh = problemData.barGlobJh + (problemData.barGlobS{n} * hJmpLambdaE0T) ./ kron(heightV0T1D(:, n-2), ones(problemData.barN, 1));
end % for n

problemData.tildeGlobP = assembleMatEdgeTrapPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDisc{1}, heightV0T1D, problemData.g.markE0Tint, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);
problemData.barGlobP = assembleMatEdge1DPhiPhiFuncDiscNuHeight(problemData.g.g1D, barU1Disc, heightV0T1D, problemData.g.g1D.markV0Tint, problemData.barHatPdiag, problemData.barHatPoffdiag);
problemData.barGlobG = assembleMatElem1DDphiPhiFuncDiscHeight(problemData.g.g1D, barU1Disc, heightQ0T1D, problemData.barHatG);

%% Assembly of boundary contributions.
u1Cont = @(x1,x2) problemData.u1Cont(t,x1,x2);
problemData.globRbdr = execin('darcyVert/assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu', problemData.g, problemData.g.markE0Tbdr, problemData.hatRdiag, DDisc);
problemData.globJD = execin('darcyVert/assembleVecEdgeTrapPhiIntFuncContNu', problemData.g, problemData.g.markE0Tbdr, u1Cont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

problemData.globPbdr = execin('darcyVert/assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu', problemData.g, problemData.g.markE0Tbdr, problemData.hatRdiag, problemData.cDisc(2:3));
problemData.globPbdr = problemData.globPbdr{1} + problemData.globPbdr{2};
problemData.barGlobPbdr = assembleMatEdge1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, heightV0T1D, problemData.g.g1D.markV0Tbdr, problemData.barHatPdiag);
end % function