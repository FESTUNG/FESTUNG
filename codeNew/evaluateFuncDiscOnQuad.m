function [dataElemOnQuad, dataEdgeIntOnQuad, dataEdgeExtOnQuad, dataEdgeExtOnQuadInt] = evaluateFuncDiscOnQuad(g, dataDisc)
global gPhi2D gPhi1D gThetaPhi1D
N = size(dataDisc, 2);
p = (sqrt(8*N+1)-3)/2; 
qOrd2D = max(2*p,1);
qOrd1D = 2*p+1;
dataEdgeIntOnQuad = cell(3,1);
dataEdgeExtOnQuad = cell(3,3);
dataEdgeExtOnQuadInt = cell(3,3);
dataElemOnQuad = dataDisc * gPhi2D{qOrd2D}.';
for nn = 1 : 3
  dataEdgeIntOnQuad{nn} = dataDisc * gPhi1D{qOrd1D}(:,:,nn).';
  for np = 1 : 3
    dataEdgeExtOnQuad{nn,np} = dataDisc * gThetaPhi1D{qOrd1D}(:,:,nn,np).';
		dataEdgeExtOnQuadInt{nn,np} = g.markE0TE0T{nn,np} * dataEdgeExtOnQuad{nn,np};
  end % for
end % for
end % function
