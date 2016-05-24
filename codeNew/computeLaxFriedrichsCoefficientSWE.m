function ret = computeLaxFriedrichsCoefficientSWE(edgeType, averaging, kronNuE0T, gConst, nn, cEdgeIntOnQuad, np, cEdgeExt2IntOnQuad, HL, uHR, vHR, HOS)
if strcmp(edgeType, 'interior')
  if strcmp(averaging, 'full-harmonic')
    assert(nargin == 9, 'Invalid number of input arguments.')
    HR = sqrt(cEdgeExt2IntOnQuad{1,nn,np});
    ret = setNaN2Zero( abs( ( ( HR .* cEdgeIntOnQuad{2,nn} + HL .* cEdgeExt2IntOnQuad{2,nn,np} ) .* kronNuE0T{nn,1}   ...
                            + ( HR .* cEdgeIntOnQuad{3,nn} + HL .* cEdgeExt2IntOnQuad{3,nn,np} ) .* kronNuE0T{nn,2} ) ...
                           ./ ( HR .* cEdgeIntOnQuad{1,nn} + HL .* cEdgeExt2IntOnQuad{1,nn,np} ) ) ) ...
          + sqrt( gConst * (cEdgeIntOnQuad{1,nn}.^1.5 + cEdgeExt2IntOnQuad{1,nn,np}.^1.5) ./ (HL + HR) );
  elseif strcmp(averaging, 'semi-harmonic')
    assert(nargin == 9, 'Invalid number of input arguments.')
    HR = sqrt(cEdgeExt2IntOnQuad{1,nn,np});
    ret = setNaN2Zero( abs( ( ( HR .* cEdgeIntOnQuad{2,nn} + HL .* cEdgeExt2IntOnQuad{2,nn,np} ) .* kronNuE0T{nn,1}   ...
                            + ( HR .* cEdgeIntOnQuad{3,nn} + HL .* cEdgeExt2IntOnQuad{3,nn,np} ) .* kronNuE0T{nn,2} ) ...
                           ./ ( HR .* cEdgeIntOnQuad{1,nn} + HL .* cEdgeExt2IntOnQuad{1,nn,np} ) ) ) ...
          + sqrt( gConst / 2 * (cEdgeIntOnQuad{1,nn} + cEdgeExt2IntOnQuad{1,nn,np}) );
  elseif strcmp(averaging, 'mean')
    assert(nargin == 8, 'Invalid number of input arguments.')
    ret = setNaN2Zero( abs( ( (cEdgeExt2IntOnQuad{1,nn,np} .* cEdgeIntOnQuad{2,nn} + cEdgeIntOnQuad{1,nn} .* cEdgeExt2IntOnQuad{2,nn,np}) .* kronNuE0T{nn,1}	 ...
                            + (cEdgeExt2IntOnQuad{1,nn,np} .* cEdgeIntOnQuad{3,nn} + cEdgeIntOnQuad{1,nn} .* cEdgeExt2IntOnQuad{3,nn,np}) .* kronNuE0T{nn,2} ) ...
                           ./ (2 * cEdgeIntOnQuad{1,nn} .* cEdgeExt2IntOnQuad{1,nn,np}) ) ) ...
          + sqrt( gConst / 2 * (cEdgeIntOnQuad{1,nn} + cEdgeExt2IntOnQuad{1,nn,np}) );
  else
    error('Unknown type of averaging.');
  end % if
elseif strcmp(edgeType, 'land')
  assert(nargin == 11, 'Invalid number of input arguments.');
  ret = ( (cEdgeIntOnQuad{2,nn} + uHR) .* kronNuE0T{nn,1} + (cEdgeIntOnQuad{3,nn} + vHR) .* kronNuE0T{nn,2} ) ./ (2*HL) + sqrt(gConst * HL);
elseif  strcmp(edgeType, 'openSea')
  assert(nargin == 12, 'Invalid number of input arguments.')
  if strcmp(averaging, 'full-harmonic')
    HR = sqrt(HOS{nn});
    ret = abs( ( cEdgeIntOnQuad{2,nn} .* kronNuE0T{nn,1} + cEdgeIntOnQuad{3,nn} .* kronNuE0T{nn,2} ) ./ (HL .* HR) ) ...
          + sqrt( gConst * (cEdgeIntOnQuad{1,nn}.^1.5 + HOS{nn}.^1.5) ./ (HL + HR) );
  elseif strcmp(averaging, 'semi-harmonic')
    ret = abs( ( cEdgeIntOnQuad{2,nn} .* kronNuE0T{nn,1} + cEdgeIntOnQuad{3,nn} .* kronNuE0T{nn,2} ) ./ (HL .* sqrt(HOS{nn})) ) ...
          + sqrt( gConst / 2 * (cEdgeIntOnQuad{1,nn} + HOS{nn}) );
  elseif strcmp(averaging, 'mean')
    ret = abs( ( cEdgeIntOnQuad{2,nn} .* kronNuE0T{nn,1} + cEdgeIntOnQuad{3,nn} .* kronNuE0T{nn,2} ) .* ( HOS{nn} + cEdgeIntOnQuad{1,nn} ) ...
            ./ ( 2 * (cEdgeIntOnQuad{1,nn} .* HOS{nn}) ) ) + sqrt( gConst / 2 * (cEdgeIntOnQuad{1,nn} + HOS{nn}) );
  else
    error('Unknown type of averaging.');
  end % if
else
  error('Riemann solver usage for any edges other then interior, land or open sea is not supported.')  
end % if
end % function
