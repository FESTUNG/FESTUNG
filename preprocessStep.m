function problemData = preprocessStep(problemData, nStep)

t = nStep * problemData.tau;

%% L2-projections of algebraic coefficients.
KDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) c(t,x1,x2), problemData.N, problemData.qOrd, ...
                       problemData.globM, problemData.basesOnQuad), problemData.KCont, 'UniformOutput', false);
fDisc = projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), problemData.N, problemData.qOrd, ...
          problemData.globM, problemData.basesOnQuad);
                             
%% Assembly of time-dependent global matrices.
problemData.globG = assembleMatElemTrapDphiPhiFuncDisc(problemData.g, problemData.hatG, KDisc);
problemData.globR = assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, ...
                    problemData.hatRdiag, problemData.hatRoffdiag, KDisc);

%% Assembly of Dirichlet boundary contributions.
hDCont = @(x1,x2) problemData.hDCont(t,x1,x2);
problemData.globRD = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, ...
                     problemData.g.markE0TbdrD, problemData.hatRdiag, KDisc);
problemData.globJD = assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                     hDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad);
problemData.globKD = problemData.eta * assembleVecEdgeTrapPhiIntFuncCont(problemData.g, ...
                     problemData.g.markE0TbdrD, hDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad);
                  
%% Assembly of Neumann boundary contributions.
gNCont = @(x1,x2) problemData.gNCont(t,x1,x2);
problemData.globKN = assembleVecEdgeTrapPhiIntFuncCont(problemData.g, ...
                     problemData.g.markE0TbdrN .* problemData.g.areaE0T, ...
                     gNCont, problemData.N, problemData.qOrd, problemData.basesOnQuad);
                   
%% Assembly of the source contribution.
problemData.globL = problemData.globM * reshape(fDisc', problemData.g.numT * problemData.N, 1);

end % function
