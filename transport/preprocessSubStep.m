function problemData = preprocessSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;
p = problemData.p;

% Project velocity to DG space, if not given
if ~all(isfield(problemData, {'hDisc', 'uHDisc', 'vHDisc', 'vNormalOnQuadEdge'})) % TODO evtl Fehler falls manche dieser Felder bekannt sind
  problemData.hDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.hCont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                2*p, problemData.hatM, problemData.basesOnQuad);
  problemData.uHDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.uHCont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                2*p, problemData.hatM, problemData.basesOnQuad);
  problemData.vHDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.vHCont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                2*p, problemData.hatM, problemData.basesOnQuad);
  % Evaluate normal velocity in quadrature points of edges
  problemData.vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.uHCont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              @(x1,x2) problemData.vHCont(problemData.timeLvls(nSubStep),x1,x2), 2*p+1);
end % if

% Evaluate cDisc in all quadrature points
qOrd2D = max(2*p,1);
problemData.cHQ0T = cellfun(@(c) (reshape(c, N, K).' * problemData.basesOnQuad.phi2D{qOrd2D}.'), problemData.cDiscRK, 'UniformOutput', false); % K x numQuad2D

% Computing the concentration from the depth-integrated one
problemData.cQ0T = cellfun(@(c) c ./ (problemData.hDisc * problemData.basesOnQuad.phi2D{qOrd2D}.'), problemData.cHQ0T, 'UniformOutput', false);
problemData.concentrationDiscRK = cellfun(@(c) problemData.swe_projectDataQ0T2DataDisc(c, 2*problemData.p, problemData.hatM, problemData.basesOnQuad), problemData.cQ0T, 'UniformOutput', false);

end % function
