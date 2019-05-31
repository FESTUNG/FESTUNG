
function pd = solveStep(pd, nStep)
KInf = pd.KInf;
NInf = pd.NInf;

% Obtain Runge-Kutta rule
switch pd.schemeType
  case 'explicit'
    %% affine vector holding all non-solution-dependent contributions
     %column:                                                                                                          HInf 1;                                                                                                             uHInf 2;                                                                                                               vHInf 3;      
     sysV = [                                                                           pd.globLInf{3} + pd.globVOSRiemInf{2};                                                           -0.5*pd.globUOSInf{1} + pd.globLInf{1} + pd.globROSInf{1};                                                            -0.5*pd.globUOSInf{2} + pd.globLInf{2} + pd.globROSInf{2}];
    %% stiffness matrix
     sysA = [                                                                              pd.globVOSRiemInf{1} + pd.globVInf,                                              pd.globQOSInf{1} + pd.globHExpInf{1} + pd.globHInf{1} + pd.globQInf{1},                                                pd.globQOSInf{2} + pd.globHExpInf{2} + pd.globHInf{2} + pd.globQInf{2};
                pd.globFInf{1} + pd.globFExpInf{1} + pd.globRInf{1} + pd.globRLInf{1} + pd.globROSRiemInf{1} + pd.globGInf{1},   pd.globOInf + pd.globOExpInf{1} + pd.globOExpInf{2} + pd.globUInf + pd.globVInf + pd.globUOSRiemInf + pd.globEInf,                                                                                                          -pd.globDInf;
                pd.globFInf{2} + pd.globFExpInf{2} + pd.globRInf{2} + pd.globRLInf{2} + pd.globROSRiemInf{2} + pd.globGInf{2},                                                                                                         pd.globDInf,   pd.globOInf + pd.globOExpInf{1} + pd.globOExpInf{2} + pd.globUInf + pd.globVInf + pd.globUOSRiemInf + pd.globEInf ];    
  

  otherwise
    error('Invalid time stepping scheme')
end % switch

pd.sysY = pd.sysY + pd.dt * (pd.sysW \ (sysV - sysA * pd.sysY));

pd.cDiscInf(:,:,1)  = reshape(pd.sysY(              1 :   KInf*NInf), NInf, KInf).';
pd.cDiscInf(:,:,2)  = reshape(pd.sysY(  KInf*NInf + 1 : 2*KInf*NInf), NInf, KInf).';
pd.cDiscInf(:,:,3)  = reshape(pd.sysY(2*KInf*NInf + 1 : 3*KInf*NInf), NInf, KInf).';
end % function