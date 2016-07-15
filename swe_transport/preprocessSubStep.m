function problemData = preprocessSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.swe_preprocessSubStep(problemData.sweData, nStep, nSubStep);

% use velocities of swe for transport problem
problemData.transportData.vNormalOnQuadEdge = problemData.sweData.massFluxQ0E0T;
problemData.transportData.u1Disc = problemData.sweData.u1Disc;
problemData.transportData.u2Disc = problemData.sweData.u2Disc;

problemData.transportData = problemData.transport_preprocessSubStep(problemData.transportData, nStep, nSubStep);
end % function
