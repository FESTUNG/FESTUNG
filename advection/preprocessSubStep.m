function problemData = preprocessSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;

% L2 projections of algebraic coefficients
fDisc  = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(problemData.t(nSubStep),x1,x2),  ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.t(nSubStep),x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.t(nSubStep),x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);

% Evaluate normal velocity in quadrature points of edges
vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(problemData.t(nSubStep),x1,x2), ...
                      @(x1,x2) problemData.u2Cont(problemData.t(nSubStep),x1,x2), 2*problemData.p+1);

% Assembly of time-dependent global matrices
problemData.globG = assembleMatElemDphiPhiFuncDiscVec(problemData.g, problemData.hatG, u1Disc, u2Disc);
problemData.globR = assembleMatEdgePhiPhiValUpwind(problemData.g, ~problemData.g.markE0TbdrN, ...
                                                   problemData.hatRdiagOnQuad, problemData.hatRoffdiagOnQuad, ...
                                                   vNormalOnQuadEdge, problemData.g.areaE0TbdrNotN);
                                               
% Assembly of Dirichlet boundary contributions
problemData.globKD = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
                        @(x1,x2) problemData.cDCont(problemData.t(nSubStep),x1,x2), vNormalOnQuadEdge, N, ...
                        problemData.basesOnQuad, problemData.g.areaE0TbdrD);

% Assembly of Neumann boundary contributions
gNUpwind = @(x1,x2) (problemData.gNCont(problemData.t(nSubStep),x1,x2) <= 0) .* problemData.gNCont(problemData.t(nSubStep),x1,x2);
problemData.globKN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                                                   gNUpwind, N, problemData.basesOnQuad);

% Assembly of the source contribution
problemData.globL = problemData.globM * reshape(fDisc', K*N, 1);
end % function
