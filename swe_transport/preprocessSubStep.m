function problemData = preprocessSubStep(problemData, nStep, nSubStep)

% TODO swe/preprocessSubStep once implemented, remove projectQuotient:
% dont't compute these twice
problemData.transportData.u1Disc = projectQuotientDisc2Disc(problemData.sweData.cDisc(:,:,2), ...
                                                            problemData.sweData.cDisc(:,:,1) - problemData.sweData.zbDisc, ...
                                                            2*problemData.sweData.p, problemData.sweData.refElemPhiPhi, problemData.sweData.basesOnQuad);
problemData.transportData.u2Disc = projectQuotientDisc2Disc(problemData.sweData.cDisc(:,:,3), ...
                                                            problemData.sweData.cDisc(:,:,1) - problemData.sweData.zbDisc, ...
                                                            2*problemData.sweData.p, problemData.sweData.refElemPhiPhi, problemData.sweData.basesOnQuad);
% use mass flux of swe in quadrature points of edges
problemData.transportData.vNormalOnQuadEdge = problemData.sweData.massFlux;
% problemData.transportData.vNormalOnQuadEdge = computeFuncDiscNuOnQuadEdge(problemData.transportData.g, problemData.transportData.u1Disc, ...
%                                                                           problemData.transportData.u2Disc, 2*problemData.transportData.p+1);

addpath('transport');
problemData.transportData = preprocessSubStep(problemData.transportData, nStep, nSubStep);
rmpath('transport');
end % function
