function [lambda, lambdaOSRiem] = computeLaxFriedrichsCoefficientsSWE(g, gConst, cEdgeInt, cEdgeExt, HOS, R1D, averaging, OSRiem)
lambda = cell(3,3);
lambdaOSRiem = cell(3,1);
for r = 1 : R1D
	for nn = 1 : 3
		if strcmp(averaging, 'full-harmonic')
			HLPow1o2 = sqrt(cEdgeInt{1,nn}(:,r));
			for np = 1 : 3
				HRPow1o2 = sqrt(cEdgeExt{1,nn,np}(:,r));
				lambda{nn,np}(:,r) = setNaN2Zero( abs( ( ( HRPow1o2 .* cEdgeInt{2,nn}(:,r) + HLPow1o2 .* cEdgeExt{2,nn,np}(:,r) ) .* g.nuE0T(:,nn,1)   ...
																							 + ( HRPow1o2 .* cEdgeInt{3,nn}(:,r) + HLPow1o2 .* cEdgeExt{3,nn,np}(:,r) ) .* g.nuE0T(:,nn,2) ) ...
																							./ ( HRPow1o2 .* cEdgeInt{1,nn}(:,r) + HLPow1o2 .* cEdgeExt{1,nn,np}(:,r) ) ) ) ...
																							+ sqrt( gConst * (cEdgeInt{1,nn}(:,r).^1.5 + cEdgeExt{1,nn,np}(:,r).^1.5) ./ (HLPow1o2 + HRPow1o2) ); % possibly NaN not necessary
% 																						find(lambda{nn,np}(:,r)== NaN)
			end % for
			if OSRiem
				HRPow1o2 = sqrt(HOS{nn}(:,r));
				lambdaOSRiem{nn}(:,r) = setNaN2Zero( abs( (cEdgeInt{2,nn}(:,r) .* g.nuE0T(:,nn,1) + cEdgeInt{3,nn}(:,r) .* g.nuE0T(:,nn,2)) ./ (HLPow1o2 .* HRPow1o2) ) ) ...
																+ sqrt( gConst * (cEdgeInt{1,nn}(:,r).^1.5 + HOS{nn}(:,r).^1.5) ./ (HLPow1o2 + HRPow1o2) );
%																find(lambdaOSRiem{nn}(:,r)== NaN)
			else
				lambdaOSRiem{nn} = [];
			end % if
		elseif strcmp(averaging, 'semi-harmonic')
			HLPow1o2 = sqrt(cEdgeInt{1,nn}(:,r));
			for np = 1 : 3
				HRPow1o2 = sqrt(cEdgeExt{1,nn,np}(:,r));
				lambda{nn,np}(:,r) = setNaN2Zero( abs( ( ( HRPow1o2 .* cEdgeInt{2,nn}(:,r) + HLPow1o2 .* cEdgeExt{2,nn,np}(:,r) ) .* g.nuE0T(:,nn,1)   ...
																							 + ( HRPow1o2 .* cEdgeInt{3,nn}(:,r) + HLPow1o2 .* cEdgeExt{3,nn,np}(:,r) ) .* g.nuE0T(:,nn,2) ) ...
																							./ ( HRPow1o2 .* cEdgeInt{1,nn}(:,r) + HLPow1o2 .* cEdgeExt{1,nn,np}(:,r) ) ) ) ...
																							+ sqrt( gConst / 2 * (cEdgeInt{1,nn}(:,r) + cEdgeExt{1,nn,np}(:,r)) ); % possibly NaN not necessary
% 																						find(lambda{nn,np}(:,r)== NaN)
			end % for
			if OSRiem
				lambdaOSRiem{nn}(:,r) = setNaN2Zero( abs( (cEdgeInt{2,nn}(:,r) .* g.nuE0T(:,nn,1) + cEdgeInt{3,nn}(:,r) .* g.nuE0T(:,nn,2)) ./ (HLPow1o2 .* sqrt(HOS{nn}(:,r))) ) ) ...
																+ sqrt( gConst / 2 * (cEdgeInt{1,nn}(:,r) + HOS{nn}(:,r)) ); % possibly NaN not necessary
%																find(lambdaOSRiem{nn}(:,r)== NaN)
			else
				lambdaOSRiem{nn} = [];
			end % if
		elseif strcmp(averaging, 'mean')
			for np = 1 :  3
				lambda{nn,np}(:,r) = setNaN2Zero( abs( ( (cEdgeExt{1,nn,np}(:,r) .* cEdgeInt{2,nn}(:,r) + cEdgeInt{1,nn}(:,r) .* cEdgeExt{2,nn,np}(:,r)) .* g.nuE0T(:,nn,1)		...
																							 + (cEdgeExt{1,nn,np}(:,r) .* cEdgeInt{3,nn}(:,r) + cEdgeInt{1,nn}(:,r) .* cEdgeExt{3,nn,np}(:,r)) .* g.nuE0T(:,nn,2) ) ...
																							./ (2 * cEdgeInt{1,nn}(:,r) .* cEdgeExt{1,nn,np}(:,r)) ) ) ...
																							+ sqrt( gConst / 2 * (cEdgeInt{1,nn}(:,r) + cEdgeExt{1,nn,np}(:,r)) ); % possibly NaN not necessary
% 																						find(lambda{nn,np}(:,r)== NaN)
			end % for
			if OSRiem
				lambdaOSRiem{nn}(:,r) = setNaN2Zero( abs( ( cEdgeInt{2,nn}(:,r) .* g.nuE0T(:,nn,1) + cEdgeInt{3,nn}(:,r) .* g.nuE0T(:,nn,2) ) ...
																								 .* (HOS{nn}(:,r) + cEdgeInt{1,nn}(:,r)) ./ (2 * (cEdgeInt{1,nn}(:,r) .* HOS{nn}(:,r))) ) ) ...
																								 + sqrt( gConst / 2 * (cEdgeInt{1,nn}(:,r) + HOS{nn}(:,r)) ); % possibly NaN not necessary
%																								 find(lambdaOSRiem{nn}(:,r)== NaN)
			else
				lambdaOSRiem{nn} = [];
			end % if
		else
			error('Unknown type of averaging.');
		end % if
	end % for
end % for
end % function
