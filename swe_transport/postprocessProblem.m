function problemData = postprocessProblem(problemData)
h = getFunctionHandle('swe/postprocessProblem');
problemData.sweData = h(problemData.sweData);

h = getFunctionHandle('transport/postprocessProblem');
problemData.transportData = h(problemData.transportData);
end % function

