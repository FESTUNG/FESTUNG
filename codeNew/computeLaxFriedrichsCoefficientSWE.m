function ret = computeLaxFriedrichsCoefficientSWE(edgeType, averaging, kronNuE0T, gConst, nn, cEdgeIntOnQuad, HEdgeIntOnQuad, np, cEdgeExt2IntOnQuad, HEdgeExt2IntOnQuad, HL, uHR, vHR, HOS)
if strcmp(edgeType, 'interior')
  if strcmp(averaging, 'full-harmonic')
    assert(nargin == 11, 'Invalid number of input arguments.')
    HR = sqrt(HEdgeExt2IntOnQuad);
    ret = setNaN2Zero( abs( ( ( HR .* cEdgeIntOnQuad{2,nn} + HL .* cEdgeExt2IntOnQuad{2,nn,np} ) .* kronNuE0T{nn,1}   ...
                            + ( HR .* cEdgeIntOnQuad{3,nn} + HL .* cEdgeExt2IntOnQuad{3,nn,np} ) .* kronNuE0T{nn,2} ) ...
                           ./ ( HR .* HEdgeIntOnQuad       + HL .* HEdgeExt2IntOnQuad ) ) ) ...
          + sqrt( gConst * (HEdgeIntOnQuad.^1.5 + HEdgeExt2IntOnQuad.^1.5) ./ (HL + HR) );
  elseif strcmp(averaging, 'semi-harmonic')
    assert(nargin == 9, 'Invalid number of input arguments.')
    HR = sqrt(HEdgeExt2IntOnQuad);
    ret = setNaN2Zero( abs( ( ( HR .* cEdgeIntOnQuad{2,nn} + HL .* cEdgeExt2IntOnQuad{2,nn,np} ) .* kronNuE0T{nn,1}   ...
                            + ( HR .* cEdgeIntOnQuad{3,nn} + HL .* cEdgeExt2IntOnQuad{3,nn,np} ) .* kronNuE0T{nn,2} ) ...
                           ./ ( HR .* HEdgeIntOnQuad       + HL .* HEdgeExt2IntOnQuad ) ) ) ...
          + sqrt( gConst / 2 * (HEdgeIntOnQuad + HEdgeExt2IntOnQuad) );
  elseif strcmp(averaging, 'mean')
    assert(nargin == 10, 'Invalid number of input arguments.')
    ret = setNaN2Zero( abs( ( (HEdgeExt2IntOnQuad .* cEdgeIntOnQuad{2,nn} + HEdgeIntOnQuad .* cEdgeExt2IntOnQuad{2,nn,np}) .* kronNuE0T{nn,1}	 ...
                            + (HEdgeExt2IntOnQuad .* cEdgeIntOnQuad{3,nn} + HEdgeIntOnQuad .* cEdgeExt2IntOnQuad{3,nn,np}) .* kronNuE0T{nn,2} ) ...
                           ./ (2 * HEdgeIntOnQuad .* HEdgeExt2IntOnQuad) ) ) ...
          + sqrt( gConst / 2 * (HEdgeIntOnQuad + HEdgeExt2IntOnQuad) );
  else
    error('Unknown type of averaging.');
  end % if
elseif strcmp(edgeType, 'land')
  assert(nargin == 13, 'Invalid number of input arguments.');
  ret = ( (cEdgeIntOnQuad{2,nn} + uHR) .* kronNuE0T{nn,1} + (cEdgeIntOnQuad{3,nn} + vHR) .* kronNuE0T{nn,2} ) ./ (2*HL) + sqrt(gConst * HL);
elseif  strcmp(edgeType, 'openSea')
  assert(nargin == 14, 'Invalid number of input arguments.')
  if strcmp(averaging, 'full-harmonic')
    HR = sqrt(HOS);
    ret = abs( ( cEdgeIntOnQuad{2,nn} .* kronNuE0T{nn,1} + cEdgeIntOnQuad{3,nn} .* kronNuE0T{nn,2} ) ./ (HL .* HR) ) ...
          + sqrt( gConst * (HEdgeIntOnQuad.^1.5 + HOS.^1.5) ./ (HL + HR) );
  elseif strcmp(averaging, 'semi-harmonic')
    ret = abs( ( cEdgeIntOnQuad{2,nn} .* kronNuE0T{nn,1} + cEdgeIntOnQuad{3,nn} .* kronNuE0T{nn,2} ) ./ (HL .* sqrt(HOS)) ) ...
          + sqrt( gConst / 2 * (HEdgeIntOnQuad + HOS) );
  elseif strcmp(averaging, 'mean')
    ret = abs( ( cEdgeIntOnQuad{2,nn} .* kronNuE0T{nn,1} + cEdgeIntOnQuad{3,nn} .* kronNuE0T{nn,2} ) .* ( HOS + HEdgeIntOnQuad ) ...
            ./ ( 2 * (HEdgeIntOnQuad .* HOS) ) ) + sqrt( gConst / 2 * (HEdgeIntOnQuad + HOS) );
  else
    error('Unknown type of averaging.');
  end % if
else
  error('Riemann solver usage for any edges other then interior, land or open sea is not supported.')  
end % if
end % function
