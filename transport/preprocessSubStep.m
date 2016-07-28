function problemData = preprocessSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;
p = problemData.p;

% Project velocity to DG space, if not given
if ~all(isfield(problemData, {'u1Disc', 'u2Disc', 'vNormalOnQuadEdge'}))
  problemData.u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                2*p, problemData.hatM, problemData.basesOnQuad);
  problemData.u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                2*p, problemData.hatM, problemData.basesOnQuad);
  % Evaluate normal velocity in quadrature points of edges
  problemData.vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              @(x1,x2) problemData.u2Cont(problemData.timeLvls(nSubStep),x1,x2), 2*p+1);
end % if

% Evaluate cDisc in all quadrature points
qOrd2D = max(2*p,1);
problemData.cQ0T = cellfun(@(c) (problemData.basesOnQuad.phi2D{qOrd2D} * reshape(c, N, K)).', problemData.cDiscRK, 'UniformOutput', false);
end % function
