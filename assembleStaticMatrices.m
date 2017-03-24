function problemData = assembleStaticMatrices(problemData)
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);

globH = assembleMatElemDphiPhi(problemData.g, problemData.hatH);
globQ = problemData.fn_assembleMatEdgeTrapPhiPhiNu(problemData.g, problemData.g.markE0Tint, problemData.hatQdiag, problemData.hatQoffdiag);
globQbdr = problemData.fn_assembleMatEdgeTrapPhiIntPhiIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU, problemData.hatQdiag);

globQAvg = problemData.fn_assembleMatEdgeTrapPhiPhiNu(problemData.g, problemData.g.markE0Tint & problemData.g.markE0Th, problemData.hatQdiag, problemData.hatQoffdiag);
globQup = assembleMatEdgeTrapPhiPhiNuBottomUp(problemData.g, problemData.g.markE0Tint | problemData.g.markE0TbdrTop, problemData.hatQdiag, problemData.hatQoffdiag);

problemData.globHQ = cellfun(@(H, Q, Qbdr) H - Q - Qbdr, globH, globQ, globQbdr, 'UniformOutput', false);
problemData.globHQup = globH{2} - globQup;
problemData.globHQavg = -globH{1} + globQAvg{1} + globQbdr{1};

problemData.globS = assembleMatEdgeTrapPhiPerQuad(problemData.g, problemData.hatSdiag);
problemData.barGlobS = assembleMatEdgeTrapPhi1DPerQuad(problemData.g, problemData.barHatSdiag);

tildeGlobH = assembleMatElemDphiPhi1D(problemData.g, problemData.tildeHatH);
tildeGlobQ = assembleMatEdgeTrapPhiPhi1DNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, problemData.tildeHatQdiag, problemData.tildeHatQoffdiag);
tildeGlobQbdr = assembleMatEdgeTrapPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrH, problemData.tildeHatQdiag);
problemData.tildeGlobHQ = problemData.gConst * (tildeGlobH{1} - tildeGlobQ{1} - tildeGlobQbdr{1});
end % function