function problemData = solveStep(problemData, ~)
K = problemData.K;
N = problemData.N;
%% Building and solving the system.
sysA = [                                                problemData.globM,                                                  sparse(K*N,K*N), -problemData.globH{1}+problemData.globQ{1}+problemData.globQN{1};
                                                          sparse(K*N,K*N),                                                problemData.globM, -problemData.globH{2}+problemData.globQ{2}+problemData.globQN{2};
         -problemData.globG{1}+problemData.globR{1}+problemData.globRD{1}, -problemData.globG{2}+problemData.globR{2}+problemData.globRD{2},                             problemData.globS+problemData.globSD];
sysV = [-problemData.globJD{1}; -problemData.globJD{2}; problemData.globKD-problemData.globKN+problemData.globL];
problemData.sysY = (problemData.sysW + problemData.tau*sysA) \ (problemData.sysW*problemData.sysY + problemData.tau*sysV);
end