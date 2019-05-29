function pd = postprocessStep(pd, nStep)
% Update time level and check for simulation end
pd.t = pd.t + pd.dt;
pd.isFinished = pd.t >= pd.tEnd;
end

