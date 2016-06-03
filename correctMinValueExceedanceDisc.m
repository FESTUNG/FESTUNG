function dataDisc = correctMinValueExceedanceDisc(dataDisc, corrSys, nStep, minTol, correctionTolerance)
if nargin < 3
	nStep = [];
end % if
if nargin < 4
	minTol = 0.001;
end % if
if nargin < 5
	correctionTolerance = Inf;
end % if

[K, N] = size(dataDisc);
p = (sqrt(8*N+1)-3)/2;
if p == 0
  dataLagr = dataDisc * phi(1,1/3,1/3);
else
  dataLagr = zeros(K,3);
  for i = 1:N
    dataLagr = dataLagr + dataDisc(:,i) * phi(i, [0, 1, 0], [0, 0, 1]);
  end % for
end % if
corr = max(minTol - dataLagr, 0);
maxCorr = max(abs(corr(:)));
if maxCorr > correctionTolerance
  [indx, indy] = find(abs(corr) == maxCorr, 1);
  error([ 'Unknown at local vertex ' num2str(indy) ' of element ' num2str(indx) ...
          ' in step ' num2str(nStep) ' is ' num2str(abs(dataLagr(indx, indy))) ...
          ' (below the minimum tolerance of ' num2str(correctionTolerance) ').' ]);
end % if
if any(corr)
  [indx, indy] = find(abs(corr) == maxCorr, 1);
  warning([ 'A maximum value of ' num2str(maxCorr) ...
            ' had to be added to unknown at local vertex ' num2str(indy) ...
            ' of element ' num2str(indx) ' in step ' num2str(nStep) '.' ]);
  if p == 0
    dataDisc = dataDisc + corr / phi(1,1/3,1/3);
  else
    dataDisc(:,1:3) = dataDisc(:,1:3) + corr / corrSys;
  end % if
end % if
end % function
