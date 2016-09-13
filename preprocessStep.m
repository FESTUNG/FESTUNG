function problemData = preprocessStep(problemData, nStep)

t = nStep * problemData.tau;

%% L2-projections of algebraic coefficients.
KDisc = cellfun(@(c) projectFuncCont2DataDiscTrap(problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
                                                  problemData.hatM{1}, problemData.basesOnQuad), ...
                problemData.KCont, 'UniformOutput', false);
fDisc = projectFuncCont2DataDiscTrap(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), problemData.qOrd, ...
                                     problemData.hatM{1}, problemData.basesOnQuad);
                             
%% Assembly of time-dependent global matrices.
problemData.globG = assembleMatElemTrapDphiPhiFuncDisc(problemData.g, problemData.hatG, KDisc);
problemData.globR = assembleMatEdgeTrapPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, ...
                                                        problemData.hatRdiag, problemData.hatRoffdiag, KDisc);

%% Assembly of Dirichlet boundary contributions.
cDCont = @(x1,x2) problemData.cDCont(t,x1,x2);
problemData.globRD = assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrD, ...
                                                                  problemData.hatRdiag, KDisc);
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
