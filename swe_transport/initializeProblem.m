function problemData = initializeProblem(problemData)
h = getFunctionHandle('swe/initializeProblem');
problemData.sweData = h(problemData.sweData);

h = getFunctionHandle('transport/initializeProblem');
problemData.transportData = h(problemData.transportData);
problemData.isFinished = problemData.sweData.isFinished || problemData.transportData.isFinished;
end % function