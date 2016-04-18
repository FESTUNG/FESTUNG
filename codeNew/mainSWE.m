function [errxi, erru, errv, erruH, errvH] = mainSWE(problem)
refinement = 0;
p          = 1;
OSRiem     = true;
scheme     = 'explicit';
fluxType	 = 'Lax-Friedrichs';
averaging	 = 'full-harmonic';
wbarOpt		 = 'off';
fileTypes  = 'vtk';

fprintf('Read user input.\n');
[ interface, name, g, zbAlg, fcAlg, gConst, NOLIBF, NWP, bottomFric, F0, F0Alg, rhsAlg, F1Alg, F2Alg, xiOSAlg, t0, tEnd, dt, numSteps, ...
	minTol, output, isVisParam, isVisu, isVisuH, isVisv, isVisvH, isVisxi, isSolAvail, xi, u, v, tidalDomain, F, fT, ramping, xiRI, uRI, vRI, ...
	rhsOSAlg, NBFR, xiOSX, xiOST, useStations, NSTAE, trianglesE, coordinatesE, NSTAV, trianglesV, coordinatesV ] = userInput(problem, refinement);
%% initial conditions
if isSolAvail % use exact solution at initial time
  xi0 = @(x1,x2  ) xi(x1,x2,t0);
   u0 = @(x1,x2  )  u(x1,x2,t0);
   v0 = @(x1,x2  )  v(x1,x2,t0);
else % cold start
  xi0 = @(x1,x2) zeros(size(x1));
   u0 = @(x1,x2) zeros(size(x1));
   v0 = @(x1,x2) zeros(size(x1));
end % if

%% globally constant parameters
K							= g.numT;                         % number of triangles
N							= nchoosek(p + 2, p);             % number of local DOFs
markE0Tint		= g.idE0T == 0;                   % [K x 3] mark local edges that are interior
markE0TbdrL		= g.idE0T == 1;										% [K x 3] mark local edges on the open sea boundary
markE0TbdrRA	= g.idE0T == 2;										% [K x 3] mark local edges on the open sea boundary
markE0TbdrRI	= g.idE0T == 3;										% [K x 3] mark local edges on the open sea boundary
markE0TbdrOS  = g.idE0T == 4;										% [K x 3] mark local edges on the open sea boundary
riverBdrs = ~isequal(markE0TbdrRI, zeros(K, 3));
g = computeDerivedGridDataSWE(g, markE0Tint, markE0TbdrL, markE0TbdrRA, markE0TbdrRI, markE0TbdrOS);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K);

[~,~,W] = quadRule2D(max(2*p,1)); R2D = length(W); % Determine the number of quadrature points on elements
[Q,W] = quadRule1D(max(2*p+1,1)); R1D = length(W); % Determine the number of quadrature points on edges
Q2X1 = @(X1,X2) g.B(:,1,1)*X1 + g.B(:,1,2)*X2 + g.coordV0T(:,1,1)*ones(size(X1));
Q2X2 = @(X1,X2) g.B(:,2,1)*X1 + g.B(:,2,2)*X2 + g.coordV0T(:,1,2)*ones(size(X1));
%% use exact approximation of the piecewise linear function zbAlg
computeBasesOnQuad(3);
zbExact = projectFuncCont2DataDisc(g, zbAlg, 2, eye(3)); % This is used for the exact linear DG representation of zb
zbEval = cell(3,1);														% This is used to calculate the total height of water from the free surface elevation for open sea boundaries
zbEvalReshaped = cell(3,1);
for n = 1:3
	[Q1, Q2] = gammaMap(n, Q);
	zbEval{n} = zbAlg(Q2X1(Q1,Q2), Q2X2(Q1,Q2));
	zbEvalReshaped{n} = reshape(zbEval{n}.', K*R1D, 1); %TODO
end % for
if riverBdrs
  xiRI = kron(xiRI, ones(R1D,1));
   uRI = kron( uRI, ones(R1D,1));
   vRI = kron( vRI, ones(R1D,1));
end % if
 if ~rhsOSAlg
  for i = 1:NBFR
    xiOSX{1,i} = xiOSX{1,i} * ones(1,R1D);
    xiOSX{2,i} = xiOSX{2,i} * ones(1,R1D);
  end % for
 end % if
%% lookup table for basis function
computeBasesOnQuad(N);
%% computation of matrices on the reference triangle
refElemPhiPhiPhi     = integrateRefElemPhiPhiPhi    (N);
refElemDphiLinPhiPhi = integrateRefElemDphiLinPhiPhi(N);
refElemDphiPhi       = integrateRefElemDphiPhi      (N);
refElemPhiPhi        = integrateRefElemPhiPhi       (N);
refEdgePhiIntPhiInt  = integrateRefEdgePhiIntPhiInt (N);
refEdgePhiIntPhiExt  = integrateRefEdgePhiIntPhiExt (N);
refElemDphiPerQuad	 = integrateRefElemDphiPerQuad	(N);
refEdgePhiIntPerQuad = integrateRefEdgePhiIntPerQuad(N);

PhiPhidiag    = integrateRefEdgePhiIntPhiIntPerQuad(N);
%% L2 projections of time-independent algebraic coefficients
fcDG = projectFuncCont2DataDisc(g, @(x1,x2) fcAlg(x1,x2), 2*p, refElemPhiPhi);
% only use zbDG for visualization and for error analysis
zbDG = projectFuncCont2DataDisc(g, @(x1,x2) zbAlg(x1,x2), 2*p, refElemPhiPhi);
if isVisParam
	zbLagrange = projectDG2Lagrange(zbExact);
	visualizeDataLagr(g, zbLagrange, 'z_b', [name, '_z_b'], 0, 'output/parameter_plotting', cd, fileTypes);
	fcLagrange = projectDG2Lagrange(fcDG);
	visualizeDataLagr(g, fcLagrange, 'f_c', [name, '_f_c'], 0, 'output/parameter_plotting', cd, fileTypes);
end % if
%% assembly of matrices corresponding to linear contributions
globD   = assembleMatElemPhiPhiFuncDisc   (g,               refElemPhiPhiPhi,     fcDG                                );
globG   = assembleMatElemPhiPhiFuncDiscLin(g,               refElemDphiLinPhiPhi, zbExact                             );
globH   = assembleMatElemDphiPhi          (g,               refElemDphiPhi                                            );
globM   = assembleMatElemPhiPhi           (g,               refElemPhiPhi                                             );
globQ   = assembleMatEdgePhiPhiNu         (g, markE0Tint,   refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, g.areaNuE0Tint  );
globQOS = assembleMatEdgePhiIntPhiIntNu   (g, markE0TbdrOS, refEdgePhiIntPhiInt,                      g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu   (g, markE0TbdrRA, refEdgePhiIntPhiInt,                      g.areaNuE0TbdrRA);
sysW = blkdiag(globM, globM, globM);

linearTerms = [   sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1},        globQ{2} + globQOS{2}  + globQRA{2} - globH{2}; 
                gConst * globG{1},															 sparse(K*N,K*N),																								 -globD;
                gConst * globG{2},																				 globD,																				 sparse(K*N,K*N) ];

clear globD globG globH globQ globQOS globQRA

%% assembly of matrices corresponding to non-linearities
% element matrices
globF		= assembleMatElemDphiPerQuad(g, refElemDphiPerQuad);
if NWP == 0
	if NOLIBF == 0
		globE = assembleMatElemPhiPhi(g, refElemPhiPhi);
	elseif NOLIBF == 1
		refElemPhiPerQuad = integrateRefElemPhiPerQuad(N);
		globE = assembleMatElemPhiPerQuad(g, refElemPhiPerQuad);
	else
		error('Invalid type of bottom friction.');
	end % if
	globE = bottomFric * globE;
elseif NWP == 1
	bottomFricDG = projectFuncCont2DataDisc(g, bottomFric, 2*p, refElemPhiPhi);
	if NOLIBF == 0
		globE = assembleMatElemPhiPhiFuncDisc(g, refElemPhiPhiPhi, bottomFricDG);
	elseif NOLIBF == 1
		refElemPhiPhiPerQuad = integrateRefElemPhiPhiPerQuad(N);
		globE = assembleMatElemPhiFuncDiscPerQuad(g, refElemPhiPhiPerQuad, bottomFricDG);
	else
		error('Invalid type of bottom friction.');
	end % if
else
	error('Invalid type of bottom friction variation.');
end % if
% Include Newtonian tidal potential in right hand sides
if tidalDomain
  rhsAlg = 0;
	for n = 1:size(F, 3)
		for i = 1:2
			for j = 1:2
				F{i,j,n} = assembleMatElemPhiPhiFuncDiscConst(g, refElemPhiPhi, F{i,j,n});
			end % for
		end % for
	end % for
end % if

% edge matrices
[globRdiag, globRoffdiag] = assembleMatEdgePhiNuPerQuad(g, markE0Tint, refEdgePhiIntPerQuad);
 globV                    = assembleMatEdgePhiPerQuad  (g,             refEdgePhiIntPerQuad);

% boundary matrices
globRL  = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrL , refEdgePhiIntPerQuad, g.areaNuE0TbdrL );
globROS = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrOS, refEdgePhiIntPerQuad, g.areaNuE0TbdrOS);
globB   = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrRA, refEdgePhiIntPerQuad, g.areaNuE0TbdrRA);
for n = 1:3
	for m = 1:2
		globRdiag{n,m} = globRdiag{n,m} + globB{n,m};
	end % for
end % for
globB = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrRI, refEdgePhiIntPerQuad, g.areaNuE0TbdrRI);

%% system matrix for correction system 
corrSys = [ phi(1,0,0) phi(1,1,0) phi(1,0,1) ; ...
            phi(2,0,0) phi(2,1,0) phi(2,0,1) ; ...
            phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

%% timestep zero
nStep =  0;
t     = t0;
rampNew = ramping(t/86400);
%% discrete data
cDG        = zeros(K,N,3);
cDG(:,:,1) = projectFuncCont2DataDisc(g, @(x1,x2)  xi0(x1,x2) - zbAlg(x1,x2)              , 2*p, refElemPhiPhi);
cDG(:,:,2) = projectFuncCont2DataDisc(g, @(x1,x2) (xi0(x1,x2) - zbAlg(x1,x2)) .* u0(x1,x2), 2*p, refElemPhiPhi);
cDG(:,:,3) = projectFuncCont2DataDisc(g, @(x1,x2) (xi0(x1,x2) - zbAlg(x1,x2)) .* v0(x1,x2), 2*p, refElemPhiPhi);
cElemOnQuad = cell(3,1);
cEdgeInt = cell(3,3);
cEdgeExt = cell(3,3,3);
cEdgeExtInt = cell(3,3,3);
%% correction for unknown c1
cDG(:,:,1) = applyMinValueExceedance2DataDisc(cDG(:,:,1), corrSys, nStep, minTol, 1000);
%% visualize initial solution
UDG = zeros(K,N,2);
if p <= 2
  if isVisu || isVisv % compute velocities
		UDG = zeros(K, N, 2);
		UDG(:,:,1) = projectQuotientDisc2DG(cDG(:,:,2), cDG(:,:,1), 2*p, refElemPhiPhi);
		UDG(:,:,2) = projectQuotientDisc2DG(cDG(:,:,3), cDG(:,:,1), 2*p, refElemPhiPhi);
		visNum = nStep / output;
		if isVisu
			uLagrange = projectDataDisc2DataLagr(UDG(:,:,1)        );
			visualizeDataLagr(g,  uLagrange,  'u_h', [name '_u' ], visNum, 'output', cd, fileTypes);
		end % if
		if isVisuH
			uHLagrange = projectDataDisc2DataLagr(cDG(:,:,2)       );
			visualizeDataLagr(g, uHLagrange, 'uH_h', [name '_uH'], visNum, 'output', cd, fileTypes);
		end % if
		if isVisv
			vLagrange = projectDataDisc2DataLagr(UDG(:,:,2)        );
			visualizeDataLagr(g,  vLagrange,  'v_h', [name '_v' ], visNum, 'output', cd, fileTypes);
		end % if
		if isVisvH
			vHLagrange = projectDataDisc2DataLagr(cDG(:,:,3)       );
			visualizeDataLagr(g, vHLagrange, 'vH_h', [name '_vH'], visNum, 'output', cd, fileTypes);
		end % if
		if isVisxi
			xiLagrange = projectDataDisc2DataLagr(cDG(:,:,1) + zbDG);
			visualizeDataLagr(g, xiLagrange, 'xi_h', [name '_xi'], visNum, 'output', cd, fileTypes);
		end % if
  end % if
end % if
%% waitbar
if strcmp(wbarOpt, 'on')
	str  = strcat( ['% done. Simulating refinement level =', ' ', num2str(refinement), ', p =', ' ', num2str(p), '.' ]);
	wbar = waitbar( 0, strcat( [ 'Time stepping:', ' ', num2str(0), str ] ) );
elseif ~strcmp(wbarOpt, 'off')
	error('Invalid option for wbarOpt. Use on for waitbar output and off to turn this feature off.')
end % if
%% save data in stations
if useStations
	elevationInStations = zeros(NSTAE, N);
	elevationValues			= zeros(numSteps+1, NSTAE);
	velocityInStations	= zeros(NSTAV, N, 2);
	velocityValues			= zeros(numSteps+1, NSTAV, 2);
end % if

%% time stepping
fprintf('\nStarting time integration from %g to %g using time step size %g (%d steps).\n\n', t0, tEnd, dt, numSteps);
tic;
while t < tEnd
  %% time update
  nStep = nStep + 1;
  t = nStep * dt;
	rampOld = rampNew;
	rampNew = ramping(t/86400);
  sysY = [ reshape(cDG(:,:,1).', K*N, 1); reshape(cDG(:,:,2).', K*N, 1); reshape(cDG(:,:,3).', K*N, 1) ];
  %% compute problem specific right hand sides and open sea boundary conditions
  if strcmp(interface, 'ADCIRC')
    [HOS, globL] = computeRightHandSidesSWE( interface, N, NBFR, xiOST, xiOSX, zbEval, g, [], t, dt, ...
																						 rhsAlg, scheme, [], [], refElemPhiPhi, globM, 0, [], gConst, ...
																						 tidalDomain, fT, F, cDG, R1D, rampOld, rampNew, riverBdrs, globB, xiRI, uRI, vRI, zbEvalReshaped);
  elseif strcmp(interface, 'none')
    [HOS, globL] = computeRightHandSidesSWE( interface, N, NBFR, xiOST, xiOSX, zbEval, g, xiOSAlg, t, dt, ...
																						 rhsAlg, scheme, F1Alg, F2Alg, refElemPhiPhi, globM, F0, F0Alg, gConst, ...
																						 tidalDomain, fT, F, cDG, R1D, rampOld, rampNew,  riverBdrs,globB, xiRI, uRI, vRI, zbEvalReshaped);
  end % if
  %% lookup table for solution
	for I = 1:3
	  [cElemOnQuad{I}, cEdgeInt(I,:), cEdgeExt(I,:,:), cEdgeExtInt(I,:,:)] = evaluateFuncDiscOnQuad(g, cDG(:,:,I)); % TODO: reshaping
	end % for
	if strcmp(fluxType, 'Lax-Friedrichs')
		[lambda, lambdaOSRiem] = computeLaxFriedrichsCoefficientsSWE(g, gConst, cEdgeInt, cEdgeExtInt, HOS, R1D, averaging, OSRiem);
	else
		error('Unknown type of flux approximation.');
	end % if
  
  %% reshaping
  for I = 1:3
		cElemOnQuad{I} = reshape(cElemOnQuad{I}.', K*R2D, 1);
		for nn = 1:3
			cEdgeInt{I,nn} = reshape(cEdgeInt{I,nn}.', K*R1D, 1);
			for np = 1:3
				cEdgeExt{I,nn,np} = reshape(cEdgeExt{I,nn,np}.', K*R1D, 1);
				cEdgeExtInt{I,nn,np} = reshape(cEdgeExtInt{I,nn,np}.', K*R1D, 1);
			end % for
		end % for
	end % for
	lambda{1,1} = reshape(lambda{1,1}.', K*R1D, 1);	lambda{1,2} = reshape(lambda{1,2}.', K*R1D, 1);	lambda{1,3} = reshape(lambda{1,3}.', K*R1D, 1);
	lambda{2,1} = reshape(lambda{2,1}.', K*R1D, 1);	lambda{2,2} = reshape(lambda{2,2}.', K*R1D, 1);	lambda{2,3} = reshape(lambda{2,3}.', K*R1D, 1);
	lambda{3,1} = reshape(lambda{3,1}.', K*R1D, 1);	lambda{3,2} = reshape(lambda{3,2}.', K*R1D, 1);	lambda{3,3} = reshape(lambda{3,3}.', K*R1D, 1);
  
  uuH = cElemOnQuad{2}.^2./cElemOnQuad{1};
  uvH = cElemOnQuad{2}.*cElemOnQuad{3}./cElemOnQuad{1};
  vvH = cElemOnQuad{3}.*cElemOnQuad{3}./cElemOnQuad{1};
  quadraticH2D = 0.5*gConst*cElemOnQuad{1}.^2;
  nonLinearity = [-globF{1}*(uuH + quadraticH2D)-globF{2}*uvH; 
                  -globF{1}*uvH-globF{2}*(vvH + quadraticH2D) ];
  
	aux = sparse(2*K*N,1);
	aux2 = sparse(2*K*N,1);
	riemann = sparse(3*K*N,1);
  for nn=1:3
    for np = 1:3
      uuH = cEdgeExt{2,nn,np}.^2 ./ cEdgeExt{1,nn,np};
			uvH = cEdgeExt{2,nn,np} .* cEdgeExt{3,nn,np} ./ cEdgeExt{1,nn,np};
			vvH = cEdgeExt{3,nn,np}.^2 ./ cEdgeExt{1,nn,np};
      quadraticH = 0.5*gConst*cEdgeExt{1,nn,np}.^2;
      nonLinearity = nonLinearity + [ globRoffdiag{nn,np,1} * (uuH + quadraticH) + globRoffdiag{nn,np,2} * uvH;
                                      globRoffdiag{nn,np,1} * uvH + globRoffdiag{nn,np,2} * (vvH + quadraticH) ];

      riemann = riemann + [ globV{nn,np} * (lambda{nn,np} .* (cEdgeInt{1,nn} - cEdgeExtInt{1,nn,np}));
														globV{nn,np} * (lambda{nn,np} .* (cEdgeInt{2,nn} - cEdgeExtInt{2,nn,np}));
														globV{nn,np} * (lambda{nn,np} .* (cEdgeInt{3,nn} - cEdgeExtInt{3,nn,np})) ];
    end % for
		uHuH = cEdgeInt{2,nn}.^2;
		uHvH = cEdgeInt{2,nn} .* cEdgeInt{3,nn};
		vHvH = cEdgeInt{3,nn}.^2;
    uuH = uHuH ./ cEdgeInt{1,nn};
		uvH = uHvH ./ cEdgeInt{1,nn};
		vvH = vHvH ./ cEdgeInt{1,nn};
    quadraticH = 0.5*gConst*cEdgeInt{1,nn}.^2;
    nonLinearity = nonLinearity + [ globRdiag{nn,1} * (uuH + quadraticH) + globRdiag{nn,2} * uvH; 
                                    globRdiag{nn,1} * uvH + globRdiag{nn,2} * (vvH + quadraticH) ];
		% land boundary contibutions
		nonLinearity = nonLinearity + [globRL{nn,1}; globRL{nn,2}] * quadraticH;
		% open sea boundary contributions
		
		if OSRiem
			nonLinearity = nonLinearity + 0.5 * [ globROS{nn,1} * (uuH + quadraticH) + globROS{nn,2} * uvH; 
																						globROS{nn,1} * uvH + globROS{nn,2} * (vvH + quadraticH) ];
			
		end % if
		
		% debug
% 		HOS{n} = ones(K,R1D);
		% debug
% 		cEdgeInt{1,nn} =  ones(K*R1D,1);
% 		cEdgeInt{2,nn} =  ones(K*R1D,1);
% 		cEdgeInt{3,nn} = zeros(K*R1D,1);
		
		HOSreshaped = reshape(HOS{n}.', K*R1D, 1);
		uuH = setNaN2Zero(cEdgeInt{2,nn}.^2 ./ HOSreshaped);
		uvH = setNaN2Zero(cEdgeInt{2,nn} .* cEdgeInt{3,nn} ./ HOSreshaped);
		vvH = setNaN2Zero(cEdgeInt{3,nn}.^2 ./ HOSreshaped);
		
		aux2 = aux2 + 0.5 * gConst * [globROS{nn,1}; globROS{nn,2}] * HOSreshaped.^2; % ROS
		aux = aux + [ globROS{nn,1} * uuH + globROS{nn,2} * uvH; globROS{nn,1} * uvH + globROS{nn,2} * vvH ]; % UOS
		
  end % for
	
% 	[(aux(aux~=0)) find(aux~=0)]
% 	find(aux~=0)
% 	length(find(aux~=0))
  
	%% compute bottom friction part
	if NOLIBF == 0
		bottomFricPart = [globE * sysY(K*N+1:2*K*N); globE * sysY(2*K*N+1:3*K*N)];
	elseif NOLIBF == 1
		bottomFricNOLIBF = (cElemOnQuad{2}.^2 +cElemOnQuad{3}.^2).^0.5 ./ cElemOnQuad{1}.^2;
		bottomFricPart = [globE * (bottomFricNOLIBF .* cElemOnQuad{2}); globE * (bottomFricNOLIBF .* cElemOnQuad{3})];
	else
		error('Invalid type of bottom friction.');
	end % if
  
  for I = 1:3
		cElemOnQuad{I} = reshape(cElemOnQuad{I}, R2D, K).';
		for nn = 1:3
			cEdgeInt{I,nn} = reshape(cEdgeInt{I,nn}, R1D, K).';
			for np = 1:3
				cEdgeExt{I,nn,np} = reshape(cEdgeExt{I,nn,np}, R1D, K).';
				cEdgeExtInt{I,nn,np} = reshape(cEdgeExtInt{I,nn,np}, R1D, K).';
			end % for
		end % for
	end % for
	lambda{1,1} = reshape(lambda{1,1}, R1D, K).';	lambda{1,2} = reshape(lambda{1,2}, R1D, K).';	lambda{1,3} = reshape(lambda{1,3}, R1D, K).';
	lambda{2,1} = reshape(lambda{2,1}, R1D, K).';	lambda{2,2} = reshape(lambda{2,2}, R1D, K).';	lambda{2,3} = reshape(lambda{2,3}, R1D, K).';
	lambda{3,1} = reshape(lambda{3,1}, R1D, K).';	lambda{3,2} = reshape(lambda{3,2}, R1D, K).';	lambda{3,3} = reshape(lambda{3,3}, R1D, K).';

	globUOSold = assembleGlobUOS(g, markE0TbdrOS, PhiPhidiag, HOS, cEdgeInt, g.areaE0TbdrOS);
% 	max(abs(blkdiag(globUOSold, globUOSold) * sysY(K*N+1:end) - aux)) %/ max( max(abs(aux)), max(abs(blkdiag(globUOSold, globUOSold) * sysY(K*N+1:end))) )
	globROSold = assembleGlobROS(g, markE0TbdrOS, HOS, N, gConst, g.areaNuE0TbdrOS);
% 	max(abs([globROSold{1}; globROSold{2}] - aux2)) / max( max(abs(aux2)), max(abs([globROSold{1}; globROSold{2}])) )
	%% building and solving the system
	switch scheme
		case 'explicit'
			if OSRiem
				%% in case of Riemann solver usage on boundaries assemble further time-dependent 
				%% global open sea boundary edge integral contributions of non-linearity, discretized explicitly
% 				globROSRiemOld = assembleGlobROSRiem(g, markE0TbdrOS, hatRdiag, cDG, gConst, g.areaNuE0TbdrOS);
% 				max(abs([globROSRiemOld{1}; globROSRiemOld{2}]*sysY(1:K*N) - aux3)) / max( max(abs(aux3)), max(abs([globROSRiemOld{1}; globROSRiemOld{2}]*sysY(1:K*N))) )
% 				globUOSRiemOld = assembleGlobUOSRiem(g, markE0TbdrOS, PhiPhidiag, cEdgeInt, g.areaE0TbdrOS);
% 				max(abs(blkdiag(globUOSRiemOld, globUOSRiemOld) * sysY(K*N+1:end) - aux3)) / max( max(abs(aux3)), max(abs(blkdiag(globUOSRiemOld, globUOSRiemOld) * sysY(K*N+1:end))) )
				globVOSRiem = assembleGlobVOSRiem(g, markE0TbdrOS, HOS, PhiPhidiag, lambdaOSRiem, g.areaE0TbdrOS);
				%% affine vector holding all non-solution-dependent contributions
				sysV = [                                     globL{1} + globVOSRiem{2};                                   globL{2} + 0.5 * globROSold{1};                           globL{3} + 0.5 * globROSold{2}     ];
				%% stiffness matrix
				sysA = [                            globVOSRiem{1},                               sparse(K*N,2*K*N) ;
									 sparse(K*N,K*N), 0.5 * globUOSold,                                                 sparse(K*N,K*N) ;
									 sparse(K*N,K*N),                                            sparse(K*N,K*N)     ,   0.5 * globUOSold  ];    
			else
				%% affine vector holding all non-solution-dependent contributions
				sysV = [                                     globL{1};                     globL{2} + globROSold{1};                     globL{3} + globROSold{2} ];
				%% stiffness matrix
				sysA = [                                        sparse(K*N,3*K*N) ;
									 sparse(K*N,K*N),   globUOSold                                 sparse(K*N,K*N) ;
									 sparse(K*N,K*N),            globRLold{2}             sparse(K*N,K*N),   globUOSold ];
			end % if
			%% calculate solution at current time via explicit Euler scheme
			sysY = sysY + dt * (sysW \ (sysV - ((sysA+linearTerms) * sysY + [sparse(K*N,1); nonLinearity + bottomFricPart] + riemann)));
		case 'semi-implicit'
			if OSRiem
				%% in case of Riemann solver usage on boundaries assemble further time-dependent 
				%% global open sea boundary edge integral contributions of non-linearity, discretized explicitly
% 				globROSRiemOld = assembleGlobROSRiem(g, markE0TbdrOS, hatRdiag, cDG, gConst, g.areaNuE0TbdrOS);
% 				globUOSRiemOld = assembleGlobUOSRiem(g, markE0TbdrOS, PhiPhidiag, cEdgeInt, g.areaE0TbdrOS);
				globVOSRiem = assembleGlobVOSRiem(g, markE0TbdrOS, HOS, PhiPhidiag, lambdaOSRiem, g.areaE0TbdrOS);
				%% affine vector holding all non-solution-dependent contributions
				sysV = [                          globL{1} + globVOSRiem{2};                           globL{2} + 0.5 * globROSold{1};                         globL{3} + 0.5 * globROSold{2}   ];
				%% stiffness matrix holding all non-linear contributions
				Aexp = [                             globVOSRiem{1},                                       sparse(K*N,K*N),                                       sparse(K*N,K*N) ;
									sparse(K*N,K*N),  0.5 * globUOSold,                                       sparse(K*N,K*N) ;
									sparse(K*N,K*N), +                              sparse(K*N,K*N),   0.5 * globUOSold ];
			else
				%% affine vector holding all non-solution-dependent contributions
				sysV = [                          globL{1};             globL{2} + globROSold{1};             globL{3} + globROSold{2}   ];
				%% stiffness matrix holding all non-linear contributions
				Aexp = [                             sparse(K*N,3*K*N)   ;
									 sparse(K*N,K*N),  globUOSold                sparse(K*N,K*N)   ;
									 sparse(K*N,K*N),            sparse(K*N,K*N),  globUOSold   ];
			end % if
			%% calculate solution at current time via semi-implicit Euler scheme
			sysY = ( sysW + dt * linearTerms ) \ ( dt * sysV - dt * (Aexp * sysY + [sparse(K*N,1); nonLinearity + bottomFricPart] + riemann) + sysW * sysY );
		otherwise
			error('Invalid scheme.')  
	end % switch
	%Euler ends
  cDG(:,:,1)  = reshape(sysY(        1 :   K*N), N, K).';
  cDG(:,:,2)  = reshape(sysY(  K*N + 1 : 2*K*N), N, K).';
  cDG(:,:,3)  = reshape(sysY(2*K*N + 1 : 3*K*N), N, K).';
  %% correction for unknown c1
	cDG(:,:,1) = applyMinValueExceedance2DataDisc(cDG(:,:,1), corrSys, nStep, minTol, 1000);
  %% visualization
  if p <= 2 && mod(nStep,output) == 0
    UDG(:,:,1) = projectQuotientDisc2DG(cDG(:,:,2), cDG(:,:,1), 2*p, refElemPhiPhi);
		UDG(:,:,2) = projectQuotientDisc2DG(cDG(:,:,3), cDG(:,:,1), 2*p, refElemPhiPhi);
		visNum = nStep / output;
		if isVisu
			uLagrange = projectDataDisc2DataLagr(UDG(:,:,1)        );
			visualizeDataLagr(g,  uLagrange,  'u_h', [name '_u' ] , visNum, 'output', cd, fileTypes);
		end % if
		if isVisuH
			uHLagrange = projectDataDisc2DataLagr(cDG(:,:,2)       );
			visualizeDataLagr(g, uHLagrange, 'uH_h', [name '_uH'], visNum, 'output', cd, fileTypes);
		end % if
		if isVisv
			vLagrange = projectDataDisc2DataLagr(UDG(:,:,2)        );
			visualizeDataLagr(g,  vLagrange,  'v_h', [name '_v' ] , visNum, 'output', cd, fileTypes);
		end % if
		if isVisvH
			vHLagrange = projectDataDisc2DataLagr(cDG(:,:,3)       );
			visualizeDataLagr(g, vHLagrange, 'vH_h', [name '_vH'], visNum, 'output', cd, fileTypes);
		end % if
		if isVisxi
			xiLagrange = projectDataDisc2DataLagr(cDG(:,:,1) + zbDG);
			visualizeDataLagr(g, xiLagrange, 'xi_h', [name '_xi'], visNum, 'output', cd, fileTypes);
		end % if
	end % if
  %% save data in stations
  if useStations
		if ~(p <= 2 && mod(nStep,output) == 0) % velocities are not available yet
			UDG(:,:,1) = projectQuotientDisc2DG(cDG(:,:,2), cDG(:,:,1), 2*p, refElemPhiPhi);
			UDG(:,:,2) = projectQuotientDisc2DG(cDG(:,:,3), cDG(:,:,1), 2*p, refElemPhiPhi);
		end % if
    uLagrange  = projectDG2Lagrange(UDG(:,:,1)       );
    vLagrange  = projectDG2Lagrange(UDG(:,:,2)       );
    xiLagrange = projectDG2Lagrange(cDG(:,:,1) + zbDG);
		for stations = 1:NSTAE
			triangles = trianglesE(stations);
			elevationInStations(stations, :) = sum(xiLagrange(triangles,:), 1) / length(triangles);
			elevationValues(nStep+1, stations) = extrapolateValue(g, triangles, coordinatesE(stations,:), elevationInStations(stations, :));
		end % for
		for stations = 1:NSTAV
			triangles = trianglesV(stations);
			velocityInStations(stations, :, 1) = sum(uLagrange(triangles,:), 1) / length(triangles);
			velocityInStations(stations, :, 2) = sum(vLagrange(triangles,:), 1) / length(triangles);
			velocityValues(nStep+1, stations, 1) = extrapolateValue(g, triangles, coordinatesV(stations,:), velocityInStations(stations, :, 1));
			velocityValues(nStep+1, stations, 2) = extrapolateValue(g, triangles, coordinatesV(stations,:), velocityInStations(stations, :, 2));
		end % for
  end % if 
  %% waitbar
	if strcmp(wbarOpt, 'on')
	  percent = round( nStep / numSteps * 100);
		wbar = waitbar( percent/100, wbar, strcat( [ 'Time stepping:', ' ', num2str(percent), str ] ) );
	end % if
end % while
fprintf(['Done. ' num2str(toc) ' seconds needed for simulation.\n']);
if strcmp(wbarOpt, 'on')
	close(wbar);
end % if
%% error computation in last timestep
if isSolAvail
	if p > 2
		% velocities are not available yet
		UDG(:,:,1) = projectQuotientDisc2DG(cDG(:,:,2), cDG(:,:,1), max(2*p,1), refElemPhiPhi);
		UDG(:,:,2) = projectQuotientDisc2DG(cDG(:,:,3), cDG(:,:,1), max(2*p,1), refElemPhiPhi);
	end % if
	errxi = computeL2Error(g, cDG(:,:,1) + zbDG, @(x1,x2) xi(x1,x2,t)															 , 2*p);
	erruH = computeL2Error(g, cDG(:,:,2)       , @(x1,x2)  u(x1,x2,t) .* (xi(x1,x2,t)-zbAlg(x1,x2)), 2*p);
	errvH = computeL2Error(g, cDG(:,:,3)       , @(x1,x2)  v(x1,x2,t) .* (xi(x1,x2,t)-zbAlg(x1,x2)), 2*p);
	erru  = computeL2Error(g, UDG(:,:,1)       , @(x1,x2)  u(x1,x2,t)															 , 2*p);
	errv  = computeL2Error(g, UDG(:,:,2)       , @(x1,x2)  v(x1,x2,t)															 , 2*p);
else
  errxi = -1; erru = -1; errv = -1; erruH = -1; errvH = -1;
end % if
if useStations
	t = t0:dt:tEnd;
	for stations = 1:NSTAE
		if stations ~= 1
			figure;
		end % if
		plot(t, elevationValues(:, stations));
		title(['free surface elevation at station ', num2str(stations)]);
	end % for
	for stations = 1:NSTAV
		figure;
		plot(t, velocityValues(:, stations, 1));
		title(['x-Component velocity at station ', num2str(stations)]);
		figure;
		plot(t, velocityValues(:, stations, 2));
		title(['y-Component velocity at station ', num2str(stations)]);
	end % for
	cd output
	save('elevationInStations', 'elevationValues');
	save( 'velocityInStations',  'velocityValues');
	cd ..
end % if
end % function
