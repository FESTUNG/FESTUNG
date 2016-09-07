function problemData = preprocessSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.swe_preprocessSubStep(problemData.sweData, nStep, nSubStep);

% use velocities of swe for transport problem
problemData.transportData.vNormalOnQuadEdge = problemData.sweData.massFluxQ0E0T;
problemData.transportData.hDisc = problemData.sweData.hDisc;
problemData.transportData.uHDisc = problemData.sweData.uHDisc;
problemData.transportData.vHDisc = problemData.sweData.vHDisc;

problemData.transportData = problemData.transport_preprocessSubStep(problemData.transportData, nStep, nSubStep);
end % function
