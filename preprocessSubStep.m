function problemData = preprocessSubStep(problemData, ~, nSubStep)
qOrd2D = max(2* problemData.p,1);

if ~isfield(problemData, 'u1Disc') && ~isfield(problemData, 'u2Disc') && ~isfield(problemData, 'vNormalOnQuadEdge')
  problemData.u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  problemData.u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  % Evaluate normal velocity in quadrature points of edges
  problemData.vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              @(x1,x2) problemData.u2Cont(problemData.timeLvls(nSubStep),x1,x2), 2*problemData.p+1);
end % if
problemData.cQ0T = cellfun(@(c) c * problemData.basesOnQuad.phi2D{qOrd2D}.', problemData.cDisc, 'UniformOutput', false);
end % function
