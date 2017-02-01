function problemData = assembleStaticMatrices(problemData)
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH = assembleMatElemDphiPhi(problemData.g, problemData.hatH);
problemData.globQ = execin('darcyVert/assembleMatEdgeTrapPhiPhiNu', problemData.g, problemData.g.markE0Tint, problemData.hatQdiag, problemData.hatQoffdiag);

globQAvg = execin('darcyVert/assembleMatEdgeTrapPhiPhiNu', problemData.g, problemData.g.markE0Tint, problemData.hatQdiag, problemData.hatQoffdiag, 1:2);
problemData.globQavg = globQAvg{1};
problemData.globQup = assembleMatEdgeTrapPhiPhiNuBottomUp(problemData.g, problemData.g.markE0Tint | problemData.g.markE0TbdrF, problemData.hatQdiag, problemData.hatQoffdiag);

problemData.globS = assembleMatEdgeTrapPhiPerQuad(problemData.g, problemData.hatSdiag);
problemData.barGlobS = assembleMatEdgeTrapPhi1DPerQuad(problemData.g, problemData.barHatSdiag);

problemData.tildeGlobH = assembleMatElemDphiPhi1D(problemData.g, problemData.tildeHatH);
problemData.tildeGlobQ = assembleMatEdgeTrapPhiPhi1DNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, problemData.tildeHatQdiag, problemData.tildeHatQoffdiag);
problemData.tildeGlobQbdr = assembleMatEdgeTrapPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tbdr, problemData.tildeHatQdiag);
for m = 1 : 2
  problemData.tildeGlobH = cellfun(@(c) problemData.gConst * c, problemData.tildeGlobH, 'UniformOutput', false);
  problemData.tildeGlobQ = cellfun(@(c) problemData.gConst * c, problemData.tildeGlobQ, 'UniformOutput', false);
  problemData.tildeGlobQbdr = cellfun(@(c) problemData.gConst * c, problemData.tildeGlobQbdr, 'UniformOutput', false);
end % for m
end

