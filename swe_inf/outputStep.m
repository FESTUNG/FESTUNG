
function pd = outputStep(pd, nStep)
%% Update waitbar.
if pd.isWaitbar
  percentDone = round( nStep / pd.numSteps * 100 );
  str  = strcat( ['% done. Simulating refinement level =', ' ', num2str(pd.refinement), ', p =', ' ', num2str(pd.p), ', pInf =', ' ', num2str(pd.pInf), '.' ]);
  pd.waitbar = waitbar(percentDone / 100, pd.waitbar, strcat( [ 'Time stepping:', ' ', num2str(percentDone), str ] ) );
end % if
end % function

