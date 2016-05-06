function dataDisc = applyMinValueExceedance2DataDisc(dataDisc, corrSys, nStep, minTol, correctionTolerance)
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
  HLagrange = dataDisc * phi(1,1/3,1/3);
else
  HLagrange = zeros(K,3);
  for i = 1:N
    HLagrange = HLagrange + dataDisc(:,i) * phi(i, [0, 1, 0], [0, 0, 1]);
  end % for
end % if
corr = max(minTol - HLagrange, 0);
normCorr = norm(corr);
if p == 0
  if ~isequal(corr, zeros(K,1))
    warning([ 'Correction of varaible c1 in vertices is needed at timestep' ' ' num2str(nStep) '. Norm of correction values:' ' ' num2str(normCorr) ] );
		if normCorr > correctionTolerance
			error('The corrections needed for variable c1 excced the maximum corrections allowed.')
		end % if
    dataDisc = dataDisc + corr / phi(1,1/3,1/3);
  end %  if
else
  if ~isequal(corr, zeros(K,3))
    warning([ 'Correction of varaible c1 in vertices is needed at timestep' ' ' num2str(nStep) '. Norm of correction values:' ' ' num2str(normCorr) ] );
    if normCorr > correctionTolerance
      error('The corrections needed for variable c1 excced the maximum corrections allowed.')
    end % if
		dataDisc(:,1:3) = dataDisc(:,1:3) + corr / corrSys;
  end % if
end % if
end % function
