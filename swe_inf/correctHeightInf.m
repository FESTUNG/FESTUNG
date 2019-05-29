function cDGInf = correctHeightInf(p, pInf, cDGInf, minTol, corrSys, nStep, beta)
[KInf, NInf] = size(cDGInf(:,:,1));
if pInf == 0
  HLagrange = cDGInf(:,:,1);% * phi(1,1/3,1/3);
else
  HLagrange = zeros(KInf,2);
  for i = 1:NInf
    HLagrange = HLagrange + cDGInf(:,i,1) * phiInf(i, p, [0, 0], [0, 1], beta);
  end % for
end % if
corr = max(minTol - HLagrange, 0);
if pInf == 0
  if ~isequal(corr, zeros(KInf,1))
    warning([ 'Correction of varaible c1 in vertices is needed at timestep' ' ' num2str(nStep) '. Norm of correction values:' ' ' num2str(norm(corr)) ] );
    cDGInf(:,:,1) = cDGInf(:,:,1) + corr;% / phi(1,1/3,1/3);
  end %  if
else
  if ~isequal(corr, zeros(KInf,2))
    warning([ 'Correction of varaible c1 in vertices is needed at timestep' ' ' num2str(nStep) '. Norm of correction values:' ' ' num2str(norm(corr)) ] );
    cDGInf(:,[1 3],1) = cDGInf(:,[1 3],1) + corr / corrSys;
  end % if
end % if
end % function
