function ret = computeLaxFriedrichsCoefficientSWE(edgeType, averaging, kronNuE0T, gConst, nn, cEdgeInt, np, cEdgeExtInt, HL, uHR, vHR, HOS)
if strcmp(edgeType, 'interior')
  if strcmp(averaging, 'full-harmonic')
    assert(nargin == 9, 'Invalid number of input arguments.')
    HR = sqrt(cEdgeExtInt{1,nn,np});
    ret = setNaN2Zero( abs( ( ( HR .* cEdgeInt{2,nn} + HL .* cEdgeExtInt{2,nn,np} ) .* kronNuE0T{nn,1}   ...
                            + ( HR .* cEdgeInt{3,nn} + HL .* cEdgeExtInt{3,nn,np} ) .* kronNuE0T{nn,2} ) ...
                           ./ ( HR .* cEdgeInt{1,nn} + HL .* cEdgeExtInt{1,nn,np} ) ) ) ...
          + sqrt( gConst * (cEdgeInt{1,nn}.^1.5 + cEdgeExtInt{1,nn,np}.^1.5) ./ (HL + HR) );
  elseif strcmp(averaging, 'semi-harmonic')
    assert(nargin == 9, 'Invalid number of input arguments.')
    HR = sqrt(cEdgeExtInt{1,nn,np});
    ret = setNaN2Zero( abs( ( ( HR .* cEdgeInt{2,nn} + HL .* cEdgeExtInt{2,nn,np} ) .* kronNuE0T{nn,1}   ...
                            + ( HR .* cEdgeInt{3,nn} + HL .* cEdgeExtInt{3,nn,np} ) .* kronNuE0T{nn,2} ) ...
                           ./ ( HR .* cEdgeInt{1,nn} + HL .* cEdgeExtInt{1,nn,np} ) ) ) ...
          + sqrt( gConst / 2 * (cEdgeInt{1,nn} + cEdgeExtInt{1,nn,np}) );
  elseif strcmp(averaging, 'mean')
    assert(nargin == 8, 'Invalid number of input arguments.')
    ret = setNaN2Zero( abs( ( (cEdgeExtInt{1,nn,np} .* cEdgeInt{2,nn} + cEdgeInt{1,nn} .* cEdgeExtInt{2,nn,np}) .* kronNuE0T{nn,1}	 ...
                            + (cEdgeExtInt{1,nn,np} .* cEdgeInt{3,nn} + cEdgeInt{1,nn} .* cEdgeExtInt{3,nn,np}) .* kronNuE0T{nn,2} ) ...
                           ./ (2 * cEdgeInt{1,nn} .* cEdgeExtInt{1,nn,np}) ) ) ...
          + sqrt( gConst / 2 * (cEdgeInt{1,nn} + cEdgeExtInt{1,nn,np}) );
  else
    error('Unknown type of averaging.');
  end % if
elseif strcmp(edgeType, 'land')
  assert(nargin == 11, 'Invalid number of input arguments.');
  ret = ( (cEdgeInt{2,nn} + uHR) .* kronNuE0T{nn,1} + (cEdgeInt{3,nn} + vHR) .* kronNuE0T{nn,2} ) ./ (2*HL) + sqrt(gConst * HL);
elseif  strcmp(edgeType, 'openSea')
  assert(nargin == 12, 'Invalid number of input arguments.')
  if strcmp(averaging, 'full-harmonic')
    HR = sqrt(HOS{nn});
    ret = abs( ( cEdgeInt{2,nn} .* kronNuE0T{nn,1} + cEdgeInt{3,nn} .* kronNuE0T{nn,2} ) ./ (HL .* HR) ) ...
          + sqrt( gConst * (cEdgeInt{1,nn}.^1.5 + HOS{nn}.^1.5) ./ (HL + HR) );
  elseif strcmp(averaging, 'semi-harmonic')
    ret = abs( ( cEdgeInt{2,nn} .* kronNuE0T{nn,1} + cEdgeInt{3,nn} .* kronNuE0T{nn,2} ) ./ (HL .* sqrt(HOS{nn})) ) ...
          + sqrt( gConst / 2 * (cEdgeInt{1,nn} + HOS{nn}) );
  elseif strcmp(averaging, 'mean')
    ret = abs( ( cEdgeInt{2,nn} .* kronNuE0T{nn,1} + cEdgeInt{3,nn} .* kronNuE0T{nn,2} ) .* ( HOS{nn} + cEdgeInt{1,nn} ) ...
            ./ ( 2 * (cEdgeInt{1,nn} .* HOS{nn}) ) ) + sqrt( gConst / 2 * (cEdgeInt{1,nn} + HOS{nn}) );
  else
    error('Unknown type of averaging.');
  end % if
else
  error('Riemann solver usage for any edges other then interior, land or open sea is not supported.')  
end % if
end % function
