function [HOS, globL] = computeRightHandSidesSWE( interface, N, NBFR, xiOST, xiOSX, zbEval, g, xiOSAlg, t, dt, ...
																									rhsAlg, scheme, F1Alg, F2Alg, hatM, globM, F0, F0Alg, gConst, ...
																									tidalDomain, fT, F, cDG, R1D, rampOld, rampNew, riverBdrs, globBRI, xiRI, uRI, vRI, zbEvalReshaped)
%% assemble source term contribution at specific time in each scheme
K = g.numT;
p = (sqrt(8*N+1)-3)/2;
globL = cell(3,1);
if rhsAlg
  switch scheme
    case 'explicit'
      F1DG = projectAlg2DG(g, @(x1,x2) F1Alg(x1,x2,t-dt), 2*p, hatM);
      F2DG = projectAlg2DG(g, @(x1,x2) F2Alg(x1,x2,t-dt), 2*p, hatM);
    case 'semi_implicit'
      F1DG = projectAlg2DG(g, @(x1,x2) F1Alg(x1,x2,t   ), 2*p, hatM);
      F2DG = projectAlg2DG(g, @(x1,x2) F2Alg(x1,x2,t   ), 2*p, hatM);
    otherwise
      error('Invalid scheme.')  
  end % switch
  globL{1} = sparse(K*N,1);
  globL{2} = globM * reshape(F1DG', K*N, 1);
  globL{3} = globM * reshape(F2DG', K*N, 1);
  if F0
    switch scheme
      case 'explicit'
        F0DG = projectAlg2DG(g, @(x1,x2) F0Alg(x1,x2,t-dt), 2*p, hatM);
      case 'semi_implicit'
        F0DG = projectAlg2DG(g, @(x1,x2) F0Alg(x1,x2,t   ), 2*p, hatM);
      otherwise
        error('Invalid scheme.')  
    end % switch
    globL{1} = globM * reshape(F0DG', K*N, 1);
  end % if
elseif tidalDomain
	globL{1} = sparse(K*N,1);
  globL{2} = sparse(K*N,K*N);
  globL{3} = sparse(K*N,K*N);
  switch scheme
		case 'explicit'
			for n = 1 : size(F,3)
        globL{2} = globL{2} + fT{1,n}(t-dt)*F{1,1,n} + fT{2,n}(t-dt)*F{1,2,n};
        globL{3} = globL{3} + fT{1,n}(t-dt)*F{2,1,n} + fT{2,n}(t-dt)*F{2,2,n};
			end % for
			H = reshape(cDG(:,:,1), K*N, 1);
			globL{2} = rampOld * globL{2} * H;
			globL{3} = rampOld * globL{3} * H;
		case 'semi_implicit'
			for n = 1 : size(F,3)
        globL{2} = globL{2} + fT{1,n}(t   )*F{1,1,n} + fT{2,n}(t   )*F{1,2,n};
        globL{3} = globL{3} + fT{1,n}(t   )*F{2,1,n} + fT{2,n}(t   )*F{2,2,n};
			end % for
			H = reshape(cDG(:,:,1), K*N, 1); % TODO schoener
			globL{2} = rampNew * globL{2} * H; % TODO
			globL{3} = rampNew * globL{3} * H;
		otherwise
			error('Invalid scheme.')
  end % switch
else
  globL{1} = sparse(K*N,1); % TODO save once before time-stepping
  globL{2} = sparse(K*N,1);
  globL{3} = sparse(K*N,1);
end % if

%% compute river boundary contributions
if riverBdrs
	HRI = cell(3,1);
	URI = cell(3,1);
	VRI = cell(3,1);
	switch scheme
		case 'explicit'
			for n = 1:3
				HRI{n} = rampOld * xiRI - zbEvalReshaped{n};
				URI{n} = rampOld * uRI;
				VRI{n} = rampOld * vRI;
			end % for
		case 'semi-implicit'
			for n = 1:3
				HRI{n} = rampNew * xiRI - zbEvalReshaped{n};
				URI{n} = rampNew * uRI;
				VRI{n} = rampNew * vRI;
			end % for
		otherwise
			error('Invalid scheme.')
	end % switch
	for n = 1:3
		globL{1} = globL{1} - globBRI{n,1} * (URI{n}.*HRI{n}) + globBRI{n,2} * (VRI{n}.*HRI{n});
		globL{2} = globL{2} - globBRI{n,1} * (URI{n}.^2.*HRI{n} + 0.5*gConst * HRI{n}.^2) + globBRI{n,2} * (URI{n}.*VRI{n}.*HRI{n});
		globL{3} = globL{3} - globBRI{n,1} * (URI{n}.*VRI{n}.*HRI{n}) + globBRI{n,2} * (VRI{n}.^2.*HRI{n} + 0.5*gConst * HRI{n}.^2);
	end % for
end % for

%% compute height on open sea boundaries
HOS = cell(3,1);
if strcmp(interface, 'ADCIRC')
	xiOS = zeros(K,R1D);
  for i = 1:NBFR
    xiOS = xiOS + xiOST{1,i}(t-dt) * xiOSX{1,i} + xiOST{2,i}(t-dt) * xiOSX{2,i};
  end % for
	% Since the open boundary condition is only used for non-linear
	% contributions we discretize it explicitly. Otherwise we would have
	% to make a distinction.
	xiOS = rampOld * xiOS;
  for n = 1:3
		HOS{n} = xiOS - zbEval{n};
  end % for
elseif strcmp(interface, 'none')
  qOrd = max(2*p, 1); [Q, ~] = quadRule1D(qOrd);
  Q2X1 = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
  Q2X2 = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
  for n = 1 : 3
    [Q1, Q2] = gammaMap(n, Q);
    HOS{n} = xiOSAlg(Q2X1(Q1, Q2), Q2X2(Q1, Q2), t-dt) - zbEval{n};
  end % for
else
  error('Invalid interface.');
end % if
end % function
