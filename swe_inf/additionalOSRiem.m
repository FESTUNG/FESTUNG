function [globROSRiem, globUOSRiem, globVOSRiem] = additionalOSRiem(g, markE0TOS, HOSAlg, hatRdiag, PhiPhidiag, cDG, gConst, t, tau)
globROSRiem = assembleGlobROSRiem(g, markE0TOS,                                 hatRdiag, cDG, gConst);
globUOSRiem = assembleGlobUOSRiem(g, markE0TOS,                               PhiPhidiag, cDG        );
globVOSRiem = assembleGlobVOSRiem(g, markE0TOS, @(x1,x2) HOSAlg(x1,x2,t-tau), PhiPhidiag, cDG, gConst);
end % function
