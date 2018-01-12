function problemData = preprocessStep(problemData, nStep)

t = nStep * problemData.tau;

%% L2-projections of algebraic coefficients.
KDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
          problemData.globM, problemData.basesOnQuad), problemData.KCont, 'UniformOutput', false);
fDisc = projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), problemData.qOrd, ...
          problemData.globM, problemData.basesOnQuad);
                             
%% Assembly of time-dependent global matrices.
problemData.globG = assembleMatElemDphiPhiFuncDisc(problemData.g, problemData.hatG, KDisc);
problemData.globR = assembleMatEdgePhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatRdiag, problemData.hatRoffdiag, KDisc);

%% Assembly of Dirichlet boundary contributions.
hDCont = @(x1,x2) problemData.hDCont(t,x1,x2);
problemData.globRD = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(problemData.g, ...
                      problemData.g.markE0TbdrD | problemData.g.markE0TbdrCoupling, problemData.hatRdiag, KDisc);
problemData.globJD = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                      hDCont, problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.globKD = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrD, ...
                      hDCont, problemData.N, problemData.basesOnQuad, problemData.qOrd, ones(problemData.g.numT, 4));
                  
%% Assembly of Neumann boundary contributions.
gNCont = @(x1,x2) problemData.gNCont(t,x1,x2);
problemData.globKN = assembleVecEdgePhiIntFuncCont(problemData.g, ...
                     problemData.g.markE0TbdrN, gNCont, problemData.N, problemData.basesOnQuad, problemData.qOrd);
                   
%% Assembly of the source contribution.
problemData.globL = problemData.globM * reshape(fDisc', problemData.g.numT * problemData.N, 1);

end % function
