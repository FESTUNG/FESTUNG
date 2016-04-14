function problemData = postprocessProblem(problemData)
K = problemData.K;
N = problemData.N;
fprintf('norm(cDisc, 1) = %g\n', norm(problemData.sysY(2*K*N+1 : 3*K*N), 1));
fprintf('norm(cDisc, 2) = %g\n', norm(problemData.sysY(2*K*N+1 : 3*K*N), 2));
fprintf('norm(cDisc, inf) = %g\n', norm(problemData.sysY(2*K*N+1 : 3*K*N), inf));
end

