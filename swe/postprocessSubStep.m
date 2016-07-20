function pd = postprocessSubStep(pd, nStep, nSubStep)
K = pd.K;
N = pd.N;

% Reshape linearized vector to solution vectors
pd.cDisc(:,:,1) = reshape(pd.cDiscRK(        1 :   K*N), N, K).';
pd.cDisc(:,:,2) = reshape(pd.cDiscRK(  K*N + 1 : 2*K*N), N, K).';
pd.cDisc(:,:,3) = reshape(pd.cDiscRK(2*K*N + 1 : 3*K*N), N, K).';

% TODO: slope limiting here

% Ensure water height doesn't fall below threshold
pd.cDisc(:,:,1) = pd.swe_correctMinValueExceedanceDisc(pd.cDisc(:,:,1), pd.sysMinValueCorrection, nStep, pd.zbLagr + pd.minTol, pd.elevTol);
pd.cDiscRK(1 : pd.K*pd.N) = reshape(pd.cDisc(:,:,1).', pd.N * pd.K, 1);

pd.isSubSteppingFinished = nSubStep >= length(pd.tLvls);
end % function
