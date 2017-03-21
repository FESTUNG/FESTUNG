function problemData = assembleStaticMatrices(problemData)
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH = assembleMatElemDphiPhi(problemData.g, problemData.hatH);
problemData.globQ = problemData.fn_assembleMatEdgeTrapPhiPhiNu(problemData.g, problemData.g.markE0Tint, problemData.hatQdiag, problemData.hatQoffdiag);
problemData.globQbdr = problemData.fn_assembleMatEdgeTrapPhiIntPhiIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU, problemData.hatQdiag);

globQAvg = problemData.fn_assembleMatEdgeTrapPhiPhiNu(problemData.g, problemData.g.markE0Tint & problemData.g.markE0Th, problemData.hatQdiag, problemData.hatQoffdiag);
problemData.globQavg = globQAvg{1};
problemData.globQup = assembleMatEdgeTrapPhiPhiNuBottomUp(problemData.g, problemData.g.markE0Tint | problemData.g.markE0TbdrTop, problemData.hatQdiag, problemData.hatQoffdiag);

problemData.globS = assembleMatEdgeTrapPhiPerQuad(problemData.g, problemData.hatSdiag);
problemData.barGlobS = assembleMatEdgeTrapPhi1DPerQuad(problemData.g, problemData.barHatSdiag);

problemData.tildeGlobH = assembleMatElemDphiPhi1D(problemData.g, problemData.tildeHatH);
problemData.tildeGlobQ = assembleMatEdgeTrapPhiPhi1DNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, problemData.tildeHatQdiag, problemData.tildeHatQoffdiag);
problemData.tildeGlobQbdr = assembleMatEdgeTrapPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrH, problemData.tildeHatQdiag);

for m = 1 : 2
  problemData.tildeGlobH{m} = problemData.gConst * problemData.tildeGlobH{m};
  problemData.tildeGlobQ{m} = problemData.gConst * problemData.tildeGlobQ{m};
  problemData.tildeGlobQbdr{m} = problemData.gConst * problemData.tildeGlobQbdr{m};
end % for m

%% Computation of bathymetry gradient.
dZbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 2);
dXbot1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 2, 1) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), 1, 1);
dXzBot = problemData.g.g1D.markT2DT * ( problemData.gConst * (dZbot1D ./ dXbot1D) );
problemData.globLzBot = kron(dXzBot, eye(problemData.N, 1));
end % function