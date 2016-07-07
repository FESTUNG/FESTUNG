function dataDisc = correctMinValueExceedanceDisc(dataDisc, corrSys, nStep, minTol, correctionTolerance)
if nargin < 3
	nStep = 0;
end % if
if nargin < 4
	minTol = 0;
end % if
if nargin < 5
	correctionTolerance = Inf;
end % if

dataV0T = projectDataDisc2DataLagr(dataDisc, 1);
corr = max(minTol - dataV0T, 0);
maxCorr = max(corr(:));

if maxCorr > correctionTolerance
  [indx, indy] = find(abs(corr) == maxCorr, 1);
  error([ 'Unknown at local vertex ' num2str(indy) ' of element ' num2str(indx) ...
          ' in step ' num2str(nStep) ' is ' num2str(dataV0T(indx, indy)) ...
          ' (below the minimum tolerance of ' num2str(correctionTolerance) ').' ]);
end % if

if any(corr(:))
  [indx, indy] = find(corr == maxCorr, 1);
  warning([ 'A maximum value of ' num2str(maxCorr) ...
            ' had to be added to unknown at local vertex ' num2str(indy) ...
            ' of element ' num2str(indx) ' in step ' num2str(nStep) '.' ]);
  if size(dataDisc,2) == 1
    dataDisc = dataDisc + max(corr, [], 2) / phi(1,1/3,1/3);
  else
    dataDisc(:,1:3) = dataDisc(:,1:3) + corr / corrSys;
  end % if
end % if
end % function
