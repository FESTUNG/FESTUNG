function [globROSRiem, globUOSRiem, globVOSRiem] = additionalOSRiemInf(gInf, HOSAlg, uOSAlg, vOSAlg, markE0TOS, hatRdiag, PhiPhidiag, cDGInf, p, gConst, t, beta, gCdiagInf, gCoffdiagInf, basesOnQuad)
globROSRiem = assembleGlobROSRiemInf(gInf,                                hatRdiag, cDGInf, gConst   );
globUOSRiem = assembleGlobUOSRiemInf(gInf,                              PhiPhidiag, cDGInf,         p, gCdiagInf);
globVOSRiem = assembleGlobVOSRiemInf(gInf, @(x1,x2) HOSAlg(x1,x2,t), @(x1,x2) uOSAlg(x1,x2,t), @(x1,x2) vOSAlg(x1,x2,t), PhiPhidiag, markE0TOS, cDGInf, gConst, p, beta, gCdiagInf, gCoffdiagInf, basesOnQuad);
end % function
