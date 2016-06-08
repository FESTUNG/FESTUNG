function checkBoundsOfDataDisc(dataDisc, lowerBound, upperBound, nStep)
dataLagr = projectDataDisc2DataLagr(dataDisc);
if min(dataLagr(:)) < lowerBound || max(dataLagr(:)) > upperBound
	error(['Unknown at step ' num2str(nStep) ' is not within the specified bounds.']); % TODO evtl location : g.V0T
end % if
end % function