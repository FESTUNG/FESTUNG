function problemData = preprocessSubStep(problemData, nStep, nSubStep)

problemData.transportData.u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                            2*problemData.p, problemData.hatM, problemData.basesOnQuad);
problemData.transportData.u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                            2*problemData.p, problemData.hatM, problemData.basesOnQuad);
% Evaluate normal velocity in quadrature points of edges
problemData.transportData.vNormalOnQuadEdge = computeFuncDiscNuOnQuadEdge(problemData.g, problemData.transportData.u1Disc, ...
                                                                          problemData.transportData.u2Disc, 2*problemData.transportData.p+1);

addpath('transport');
problemData.transportData = preprocessSubStep(problemData.transportData, nStep, nSubStep);
rmpath('transport');
end % function
