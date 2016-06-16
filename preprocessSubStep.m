function problemData = preprocessSubStep(problemData, ~, nSubStep)

if ~isfield(problemData, 'u1Disc') && ~isfield(problemData, 'u2Disc') && ~isfield(problemData, 'vNormalOnQuadEdge')
  problemData.u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  problemData.u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  % Evaluate normal velocity in quadrature points of edges
  problemData.vNormalOnQuadEdge = computeFuncDiscNuOnQuadEdge(problemData.g, problemData.u1Disc, ...
                                                              problemData.u2Disc, 2*problemData.p+1);
end % if
end % function
