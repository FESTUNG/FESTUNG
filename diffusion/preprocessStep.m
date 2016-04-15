function problemData = preprocessStep(problemData, nStep)
t = nStep * problemData.tau;
%% L2-projections of algebraic coefficients.
dDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.dCont(t,x1,x2), 2 * problemData.p, ...
                                 problemData.hatM, problemData.basesOnQuad);
fDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), 2 * problemData.p, ...
                                 problemData.hatM, problemData.basesOnQuad);
%% Assembly of time-dependent global matrices.
problemData.globG = assembleMatElemDphiPhiFuncDisc(problemData.g, problemData.hatG, dDisc);
problemData.globR = assembleMatEdgePhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatRdiag, problemData.hatRoffdiag, dDisc, problemData.g.areaNuE0Tint);
%% Assembly of Dirichlet boundary contributions.
cDCont = @(x1,x2) problemData.cDCont(t,x1,x2);
problemData.globRD = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrD, ...
                      problemData.hatRdiag, dDisc, problemData.g.areaNuE0TbdrD);
problemData.globJD = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                      cDCont, problemData.N, problemData.basesOnQuad, problemData.g.areaNuE0TbdrD);
problemData.globKD = problemData.eta * assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrD, ...
                      cDCont, problemData.N, problemData.basesOnQuad);
%% Assembly of Neumann boundary contributions.
problemData.globKN = assembleVecEdgePhiIntFuncDiscIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                      dDisc, @(x1,x2) problemData.gNCont(t,x1,x2), problemData.basesOnQuad, problemData.g.areaE0TbdrN);
%% Assembly of the source contribution.
problemData.globL = problemData.globM * reshape(fDisc', problemData.K * problemData.N, 1);
end % function