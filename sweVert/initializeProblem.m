function problemData = initializeProblem(problemData)
problemData.isFinished = false;
K = problemData.g.numT;
N = problemData.N;

% TODO

% Vector of unknowns (H, U, W)
problemData.cDisc = { zeros(K,N); zeros(K,N); zeros(K,N) };
end % function