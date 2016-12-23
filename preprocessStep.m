function problemData = preprocessStep(problemData, nStep)
t = nStep * problemData.tau;

%% L2-projections of algebraic coefficients.
DDisc = cellfun(@(c) execin('darcyVert/projectFuncCont2DataDiscTrap', problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
        problemData.hatM{1}, problemData.basesOnQuad), problemData.DCont, 'UniformOutput', false);

% TODO

%% Assembly of time-dependent global matrices.
problemData.globE = execin('darcyVert/assembleMatElemTrapDphiPhiFuncDisc', problemData.g, problemData.hatG, problemData.cDisc(2:3));
problemData.globG = execin('darcyVert/assembleMatElemTrapDphiPhiFuncDisc', problemData.g, problemData.hatG, DDisc);

problemData.globR = execin('darcyVert/assembleMatEdgeTrapPhiPhiFuncDiscNu', problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);
end % function