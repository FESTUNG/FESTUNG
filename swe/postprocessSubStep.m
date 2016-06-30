function pd = postprocessSubStep(pd, nStep, nSubStep)
% Ensure water height doesn't fall below threshold
xi = reshape(pd.cDiscRK(1 : pd.K*pd.N), pd.N, pd.K).';
xi = correctMinValueExceedanceDisc(xi, pd.sysMinValueCorrection, nStep, pd.zbLagr + pd.minTol, 20);
pd.cDiscRK(1 : pd.K*pd.N) = reshape(xi.', pd.N * pd.K, 1);

pd.isSubSteppingFinished = nSubStep >= length(pd.tLvls);
end % function
