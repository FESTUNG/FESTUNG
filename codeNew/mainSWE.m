function [errxi, erru, errv, erruH, errvH] = mainSWE(problem)
global gPhi2D gPhi1D gThetaPhi1D
refinement = 0;
p          = 1;
OSRiem     = 1;
land       = 'riemann';
scheme     = 'explicit';
fluxType	 = 'Lax-Friedrichs';
averaging	 = 'full-harmonic';
wbarOpt		 = 'off';
fileTypes  = 'vtk';

fprintf('Read user input.\n');
[ interface, name, g, zbAlg, fcAlg, gConst, NOLIBF, NWP, bottomFric, F0, F0Alg, rhsAlg, F1Alg, F2Alg, xiOSAlg, t0, tEnd, dt, numSteps, ...
	minTol, output, isVisParam, isVisu, isVisuH, isVisv, isVisvH, isVisxi, isSolAvail, xi, u, v, tidalDomain, F, fT, NRAMP, ramping, ...
  xiRI, uRI, vRI, rhsOSAlg, NBFR, xiOSX, xiOST, useStations, NSTAE, triE, coordE, NSTAV, triV, coordV ] = userInput(problem, refinement);
%% initial conditions
if isSolAvail % use exact solution at initial time
  xi0 = @(x1,x2) xi(x1,x2,t0);
   u0 = @(x1,x2)  u(x1,x2,t0);
   v0 = @(x1,x2)  v(x1,x2,t0);
else % cold start
  xi0 = @(x1,x2) zeros(size(x1));
  if problem == 4
   u0 = @(x1,x2) x1==x1;
  else
    u0 = @(x1,x2) zeros(size(x1));
  end % if
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
g             = computeDerivedGridDataSWE(g, markE0Tint, markE0TbdrL, markE0TbdrRA, markE0TbdrRI, markE0TbdrOS);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K);

landBdrs      = ~isequal(markE0TbdrL , zeros(K, 3));
radiationBdrs = ~isequal(markE0TbdrRA, zeros(K, 3));
riverBdrs     = ~isequal(markE0TbdrRI, zeros(K, 3));
openSeaBdrs   = ~isequal(markE0TbdrOS, zeros(K, 3));

[~,~,W] = quadRule2D(max(2*p,1)); R2D = length(W); % Determine the number of quadrature points on elements
[Q,W] = quadRule1D(max(2*p+1,1)); R1D = length(W); % Determine the number of quadrature points on edges

if strcmp(land, 'reflected')
  nuE0Tprod = cell(3,1);
  nuE0Tsqrd = cell(3,2);
  for nn = 1:3
    nuE0Tprod{nn}   = kron(g.nuE0T(:,nn,1) .* g.nuE0T(:,nn,2), ones(R1D,1));
    nuE0Tsqrd{nn,1} = kron(g.nuE0T(:,nn,1).^2, ones(R1D,1));
    nuE0Tsqrd{nn,2} = kron(g.nuE0T(:,nn,2).^2, ones(R1D,1));
  end % for
elseif strcmp(land, 'riemann')
  nuE0Tdiff = cell(3,1);
  nuE0Tprod = cell(3,1);
  for nn = 1:3
    nuE0Tdiff{nn}   = kron(g.nuE0T(:,nn,2).^2 - g.nuE0T(:,nn,1).^2, ones(R1D,1));
    nuE0Tprod{nn}   = kron(g.nuE0T(:,nn,1) .* g.nuE0T(:,nn,2), ones(R1D,1));
  end % for
elseif ~strcmp(land, 'natural')
  error('Unknown type of land boundary discretization.')
end % if
if ~openSeaBdrs
  OSRiem = false;
end % if

%% L2 projections of time-independent algebraic coefficients
% use exact approximation of the piecewise linear function zbAlg
computeBasesOnQuad(3);
refElemPhiLinPhiLin = integrateRefElemPhiPhi(3);
fcDG = projectFuncCont2DataDisc(g, @(x1,x2) fcAlg(x1,x2), 2, refElemPhiLinPhiLin);
zbExact = projectFuncCont2DataDisc(g, zbAlg, 2, refElemPhiLinPhiLin); % This is used for the exact linear DG representation of zb
zbEvalOnQuad1D = cell(3,1);                                           % This is used for boundary conditions

if isVisParam
	zbLagrange = projectDataDisc2DataLagr(zbExact);
	visualizeDataLagr(g, zbLagrange, 'z_b', [name, '_z_b'], [], 'output', cd, fileTypes);
	fcLagrange = projectDataDisc2DataLagr(fcDG);
	visualizeDataLagr(g, fcLagrange, 'f_c', [name, '_f_c'], [], 'output', cd, fileTypes);
end % if

quadPhysPts1D = cell(3,2);
kronNuE0T = cell(3,2);
for nn = 1:3
	[Q1, Q2] = gammaMap(nn, Q);
  aux = cell(2,1);
  for m = 1:2
    aux{m} = g.mapRef2Phy(m, Q1, Q2);
    quadPhysPts1D{nn,1} = reshape(aux{m}.', R1D*K, 1);
    kronNuE0T{nn,m} = kron(g.nuE0T(:,nn,m), ones(R1D,1));
  end % for
	zbEvalOnQuad1D{nn} = zbAlg(aux{1}, aux{2});
  zbEvalOnQuad1D{nn} = reshape(zbEvalOnQuad1D{nn}.', R1D*K, 1);
end % for
if riverBdrs
  xiRI = kron(xiRI, ones(R1D,1));
   uRI = kron( uRI, ones(R1D,1));
   vRI = kron( vRI, ones(R1D,1));
end % if
if openSeaBdrs && ~rhsOSAlg
  for i = 1:NBFR
    xiOSX{1,i} = kron(xiOSX{1,i}, ones(R1D,1));
    xiOSX{2,i} = kron(xiOSX{2,i}, ones(R1D,1));
  end % for
end % if
 
%% lookup table for basis function
computeBasesOnQuad(N);
%% computation of matrices on the reference triangle
refEdgePhiIntPerQuad = integrateRefEdgePhiIntPerQuad(N);
refEdgePhiIntPhiExt  = integrateRefEdgePhiIntPhiExt (N);
refEdgePhiIntPhiInt  = integrateRefEdgePhiIntPhiInt (N);
refElemDphiPerQuad	 = integrateRefElemDphiPerQuad	(N);
refElemDphiPhi       = integrateRefElemDphiPhi      (N);
refElemPhiPhi        = integrateRefElemPhiPhi       (N);
refElemPhiPhiDphiLin = integrateRefElemPhiPhiDphiLin(N);
refElemPhiPhiPhi     = integrateRefElemPhiPhiPhi    (N);
refElemPhiPhiPhiLin  = integrateRefElemPhiPhiPhiLin (N);

% only use zbDG for visualization and for error analysis
zbDG = projectFuncCont2DataDisc(g, @(x1,x2) zbAlg(x1,x2), 2*p, refElemPhiPhi);

%% assembly of matrices corresponding to linear contributions
globQOS = assembleMatEdgePhiIntPhiIntNu    (g, markE0TbdrOS, refEdgePhiIntPhiInt,                      g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu    (g, markE0TbdrRA, refEdgePhiIntPhiInt,                      g.areaNuE0TbdrRA);
globQ   = assembleMatEdgePhiPhiNu          (g, markE0Tint,   refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, g.areaNuE0Tint  );
globH   = assembleMatElemDphiPhi           (g,               refElemDphiPhi                                            );
globM   = assembleMatElemPhiPhi						 (g,               refElemPhiPhi																						 );
globG   = assembleMatElemPhiPhiDfuncDiscLin(g,               refElemPhiPhiDphiLin, zbExact                             );
globD   = assembleMatElemPhiPhiFuncDiscLin (g,               refElemPhiPhiPhiLin,  fcDG                                );
sysW = blkdiag(globM, globM, globM);

linTerms = [   sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1}, globQ{2} + globQOS{2}  + globQRA{2} - globH{2}; 
             gConst * globG{1},															 sparse(K*N,K*N),																					 -globD;
             gConst * globG{2},																				 globD,                                  sparse(K*N,K*N) ];

clear globD globG globH globQ globQOS globQRA

%% assembly of matrices corresponding to non-linearities
% element matrices
% bottom friction matrix
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

% nonlinearities
globF		= assembleMatElemDphiPerQuad(g, refElemDphiPerQuad);

% Newtonian tide potential
if tidalDomain
  rhsAlg = 0;
	for l = 1:size(F, 3)
		for i = 1:2
			for j = 1:2
				F{i,j,l} = assembleMatElemPhiPhiFuncDiscConst(g, refElemPhiPhi, F{i,j,l});
			end % for
		end % for
	end % for
end % if

% edge matrices
[globRdiag, globRoffdiag] = assembleMatEdgePhiNuPerQuad(g, markE0Tint, refEdgePhiIntPerQuad);
 globV                    = assembleMatEdgePhiPerQuad  (g,             refEdgePhiIntPerQuad);

% boundary matrices
if landBdrs
  globRL  = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrL , refEdgePhiIntPerQuad, g.areaNuE0TbdrL);
  if strcmp(land, 'riemann')
    globVL = assembleMatEdgePhiIntPerQuad(g, markE0TbdrL, refEdgePhiIntPerQuad, g.areaE0TbdrL);
  end % if
end % if
if radiationBdrs
  globRRA   = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrRA, refEdgePhiIntPerQuad, g.areaNuE0TbdrRA);
  for nn = 1:3
    for m = 1:2
      globRdiag{nn,m} = globRdiag{nn,m} + globRRA{nn,m};
    end % for
  end % for
  clear globRRA
end % if
if riverBdrs
  globRRI = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrRI, refEdgePhiIntPerQuad, g.areaNuE0TbdrRI);
end % if
if openSeaBdrs
  HOS = cell(3,1);
  globROS = assembleMatEdgePhiIntNuPerQuad(g, markE0TbdrOS, refEdgePhiIntPerQuad, g.areaNuE0TbdrOS);
  if OSRiem
    globVOS = assembleMatEdgePhiIntPerQuad(g, markE0TbdrOS, refEdgePhiIntPerQuad, g.areaE0TbdrOS);
  end % if
end % if

globLRI = cell(3,1);
globLRI{1} = sparse(K*N,1);
globLRI{2} = sparse(K*N,1);
globLRI{3} = sparse(K*N,1);
if riverBdrs && ~NRAMP
  for nn = 1:3
      HRI = xiRI - zbEvalOnQuad1D{nn};
     uHRI = uRI .* HRI;
     vHRI = vRI .* HRI;
    uvHRI = uRI .* vHRI;
    quadraticH = 0.5 * gConst * HRI.^2;
    globLRI{1} = globLRI{1} + globRRI{nn,1} * uHRI + globRRI{nn,2} * vHRI;
    globLRI{2} = globLRI{2} + globRRI{nn,1} * (uRI .* uHRI + quadraticH) + globRRI{nn,2} * uvHRI;
    globLRI{3} = globLRI{3} + globRRI{nn,1} * uvHRI + globRRI{nn,2} * (vRI .* vHRI + quadraticH);
  end % for
end % if

% cell for algebraic right hand sides
globL = cell(3,1);
globL{1} = sparse(K*N,1);
globL{2} = sparse(K*N,1);
globL{3} = sparse(K*N,1);
% cell for tidal potential
tidePot = cell(2,1);
tidePot{1} = sparse(K*N,K*N);
tidePot{2} = sparse(K*N,K*N);

%% system matrix for correction system 
corrSys = [ phi(1,0,0) phi(1,1,0) phi(1,0,1); ...
            phi(2,0,0) phi(2,1,0) phi(2,0,1); ...
            phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

%% timestep zero
nStep =  0;
t     = t0;
if NRAMP
  rampNew = ramping(t/86400);
  rampOld = rampNew;
else
  rampNew = [];
  rampOld = [];
end % if
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
if isVisu || isVisv
  % compute velocities
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
  % time update
  nStep = nStep + 1;
  tOld  = t;
  t     = nStep * dt;
  if NRAMP
  	rampOld = rampNew;
    rampNew = ramping(t/86400);
  end % if
  % solution at old time step
  sysY = [ reshape(cDG(:,:,1).', K*N, 1); reshape(cDG(:,:,2).', K*N, 1); reshape(cDG(:,:,3).', K*N, 1) ];
  
  %% right hand sides / source terms / non solution dependent contributions
  if tidalDomain
    tidePot{1} = sparse(K*N,K*N);
		tidePot{2} = sparse(K*N,K*N);
		switch scheme
			case 'explicit'
				for n = 1 : size(F,3)
					tidePot{1} = tidePot{1} + fT{1,n}(t-dt)*F{1,1,n} + fT{2,n}(t-dt)*F{1,2,n};
					tidePot{2} = tidePot{2} + fT{1,n}(t-dt)*F{2,1,n} + fT{2,n}(t-dt)*F{2,2,n};
				end % for
				if NRAMP
					tidePot{1} = rampOld * tidePot{1};
					tidePot{2} = rampOld * tidePot{2};
				end % if
			case 'semi_implicit'
				for n = 1 : size(F,3)
					tidePot{1} = tidePot{1} + fT{1,n}(t   )*F{1,1,n} + fT{2,n}(t   )*F{1,2,n};
					tidePot{2} = tidePot{2} + fT{1,n}(t   )*F{2,1,n} + fT{2,n}(t   )*F{2,2,n};
				end % for
				if NRAMP
					tidePot{1} = rampNew * tidePot{1};
					tidePot{2} = rampNew * tidePot{2};
				end % if
			otherwise
				error('Invalid scheme.')
		end % switch
  end % if
  if rhsAlg
    if F0
      switch scheme
        case 'explicit'
          F0DG = projectFuncCont2DataDisc(g, @(x1,x2) F0Alg(x1,x2,tOld), 2*p, refElemPhiPhi);
        case 'semi_implicit'
          F0DG = projectFuncCont2DataDisc(g, @(x1,x2) F0Alg(x1,x2,t   ), 2*p, refElemPhiPhi);
        otherwise
          error('Invalid scheme.')
      end % switch
      globL{1} = globM * reshape(F0DG', K*N, 1);
    end % if
    switch scheme
      case 'explicit'
        F1DG = projectFuncCont2DataDisc(g, @(x1,x2) F1Alg(x1,x2,tOld), 2*p, refElemPhiPhi);
        F2DG = projectFuncCont2DataDisc(g, @(x1,x2) F2Alg(x1,x2,tOld), 2*p, refElemPhiPhi);
      case 'semi_implicit'
        F1DG = projectFuncCont2DataDisc(g, @(x1,x2) F1Alg(x1,x2,t   ), 2*p, refElemPhiPhi);
        F2DG = projectFuncCont2DataDisc(g, @(x1,x2) F2Alg(x1,x2,t   ), 2*p, refElemPhiPhi);
      otherwise
        error('Invalid scheme.')  
    end % switch
    globL{2} = globM * reshape(F1DG', K*N, 1);
    globL{3} = globM * reshape(F2DG', K*N, 1);
  end % if
  
  if riverBdrs && NRAMP
    switch scheme
      case 'explicit'
        xi = rampOld * xiRI;
         u = rampOld *  uRI;
         v = rampOld *  vRI;
      case 'semi-implicit'
        xi = rampNew * xiRI;
         u = rampNew *  uRI;
         v = rampNew *  vRI;
      otherwise
        error('Invalid scheme.')
    end % switch
  end % if
  
  if openSeaBdrs
    % Since the open boundary condition is only used for non-linear
    % contributions we discretize it explicitly. Otherwise we would have
    % to make a distinction.
    if strcmp(interface, 'ADCIRC')
      xiOS = zeros(K*R1D, 1);
      for i = 1 : size(xiOSX,2)
        xiOS = xiOS + xiOST{1,i}(tOld) * xiOSX{1,i} + xiOST{2,i}(tOld) * xiOSX{2,i};
      end % for
      if NRAMP
        xiOS = rampOld * xiOS;
      end % if
    end % if
  end % if
  
  % initialize fields
  riem   = sparse(3*K*N,1);
  riemOS = sparse(K*N,1);
  if riverBdrs && NRAMP
    globLRI{1} = sparse(K*N,1);
    globLRI{2} = sparse(K*N,1);
    globLRI{3} = sparse(K*N,1);
  end % if
  
  %% solution dependent contributions
  % lookup table for solution
	for I = 1:3
    cElemOnQuad{I} = reshape(gPhi2D{2*p} * cDG(:,:,I).', R2D*K, 1);
	end % for
	
  %% element contributions
  % bottom friction part
	if NOLIBF == 0
		bottomFricPart = [globE * sysY(K*N+1:2*K*N); globE * sysY(2*K*N+1:3*K*N)];
	elseif NOLIBF == 1
		bottomFricNOLIBF = (cElemOnQuad{2}.^2 +cElemOnQuad{3}.^2).^0.5 ./ cElemOnQuad{1}.^2;
		bottomFricPart = [globE * (bottomFricNOLIBF .* cElemOnQuad{2}); globE * (bottomFricNOLIBF .* cElemOnQuad{3})];
	else
		error('Invalid type of bottom friction.');
	end % if
  
  % nonlinearity
  uuH = cElemOnQuad{2}.^2 ./ cElemOnQuad{1};
  uvH = cElemOnQuad{2}.*cElemOnQuad{3} ./ cElemOnQuad{1};
  vvH = cElemOnQuad{3}.*cElemOnQuad{3} ./ cElemOnQuad{1};
  quadraticH = 0.5 * gConst * cElemOnQuad{1}.^2;
  nonLinearity = [-globF{1} * (uuH + quadraticH)-globF{2} * uvH; 
                  -globF{1} * uvH-globF{2} * (vvH + quadraticH) ];

  % edge contributions
  for nn=1:3
    for I = 1:3
      cEdgeInt{I,nn} = reshape(gPhi1D{2*p+1}(:,:,nn) * cDG(:,:,I).', R1D*K, 1);
    end % for
		if strcmp(averaging, 'full-harmonic') || strcmp(averaging, 'semi-harmonic')
			HL = sqrt(cEdgeInt{1,nn});
		end % if
		for np = 1:3
      for I = 1:3
        aux = gThetaPhi1D{2*p+1}(:,:,nn,np) * cDG(:,:,I).';
        cEdgeExt{I,nn,np} = reshape(aux, R1D*K, 1);
        cEdgeExtInt{I,nn,np} = reshape(aux * g.markE0TE0T{nn,np}.', R1D*K, 1);
      end % for
      uuH = cEdgeExt{2,nn,np}.^2 ./ cEdgeExt{1,nn,np};
			uvH = cEdgeExt{2,nn,np} .* cEdgeExt{3,nn,np} ./ cEdgeExt{1,nn,np};
			vvH = cEdgeExt{3,nn,np}.^2 ./ cEdgeExt{1,nn,np};
      quadraticH = 0.5*gConst*cEdgeExt{1,nn,np}.^2;
      nonLinearity = nonLinearity + [ globRoffdiag{nn,np,1} * (uuH + quadraticH) + globRoffdiag{nn,np,2} * uvH;
                                      globRoffdiag{nn,np,1} * uvH + globRoffdiag{nn,np,2} * (vvH + quadraticH) ];
			if strcmp(fluxType, 'Lax-Friedrichs')
        lambda = computeLaxFriedrichsCoefficientSWE( 'interior', averaging, nn, np, kronNuE0T, [], [], [], [], [], HL, [], ...
                                                     cEdgeInt, cEdgeExtInt, gConst );
				riem = riem + [ globV{nn,np} * (lambda .* (cEdgeInt{1,nn} - cEdgeExtInt{1,nn,np}));
												globV{nn,np} * (lambda .* (cEdgeInt{2,nn} - cEdgeExtInt{2,nn,np}));
												globV{nn,np} * (lambda .* (cEdgeInt{3,nn} - cEdgeExtInt{3,nn,np})) ];
			else
				error('Unknown type of flux approximation.');
			end % if
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
    if strcmp(land, 'natural')
      nonLinearity = nonLinearity + [globRL{nn,1}; globRL{nn,2}] * quadraticH;
    elseif strcmp(land, 'reflected')
      uH = nuE0Tsqrd{nn,2} .* cEdgeInt{2,nn} - nuE0Tprod{nn} .* cEdgeInt{3,nn};
      vH = nuE0Tsqrd{nn,1} .* cEdgeInt{1,nn} - nuE0Tprod{nn} .* cEdgeInt{2,nn};
      uuH = uH.^2 ./ cEdgeInt{1,nn};
      uvH = uH .* vH ./ cEdgeInt{1,nn};
      vvH = vH.^2 ./ cEdgeInt{1,nn};
      nonLinearity = nonLinearity + [ globRL{nn,1} * (uuH + quadraticH) + globRL{nn,2} * uvH;
                                      globRL{nn,1} * uvH + globRL{nn,2} * (vvH + quadraticH) ];
    elseif strcmp(land, 'riemann')
      uHR =  nuE0Tdiff{nn} .* cEdgeInt{2,nn} - 2 * nuE0Tprod{nn} .* cEdgeInt{3,nn};
      vHR = -nuE0Tdiff{nn} .* cEdgeInt{3,nn} - 2 * nuE0Tprod{nn} .* cEdgeInt{2,nn};
      uuHR = uHR.^2 ./ cEdgeInt{1,nn};
      uvHR = uHR .* vHR ./ cEdgeInt{1,nn};
      vvHR = vHR.^2 ./ cEdgeInt{1,nn};
			if strcmp(fluxType, 'Lax-Friedrichs')
	      lambda = computeLaxFriedrichsCoefficientSWE( 'land', [], nn, [], kronNuE0T, cEdgeInt{2,nn}, uHR, cEdgeInt{3,nn}, vHR, cEdgeInt{1,nn}, ...
                                                     [], [], [], [], gConst );
				nonLinearity = nonLinearity + 0.5 * ...
       [ globRL{nn,1} * (uuH + uuHR + 2 * quadraticH) + globRL{nn,2} * (uvH + uvHR) + globVL{nn} * (lambda .* (cEdgeInt{2,nn} - uHR)); ...
         globRL{nn,1} * (uvH + uvHR) + globRL{nn,2} * (vvH + vvHR + 2 * quadraticH) + globVL{nn} * (lambda .* (cEdgeInt{3,nn} - vHR)) ];
			else
				error('Unknown type of flux approximation.');
			end % if
    else
      error('Unknown type of land boundary discretization.')
    end % if
		
    % river boundary contributions
    if riverBdrs && NRAMP
        HRI = xi - zbEvalOnQuad1D{nn};
       uHRI = u .* HRI;
       vHRI = v .* HRI;
      uvHRI = u .* vHRI;
      quadraticH = 0.5 * gConst * HRI.^2;
      globLRI{1} = globLRI{1} + globRRI{nn,1} * uHRI + globRRI{nn,2} * vHRI;
      globLRI{2} = globLRI{2} + globRRI{nn,1} * (uRI .* uHRI + quadraticH) + globRRI{nn,2} * uvHRI;
      globLRI{3} = globLRI{3} + globRRI{nn,1} * uvHRI + globRRI{nn,2} * (vRI .* vHRI + quadraticH);
    end % if
    
		% open sea boundary contributions
    if OSRiem
			nonLinearity = nonLinearity + 0.5 * [ globROS{nn,1} * (uuH + quadraticH) + globROS{nn,2} * uvH; 
																						globROS{nn,1} * uvH + globROS{nn,2} * (vvH + quadraticH) ];
      if strcmp(interface, 'ADCIRC')
        HOS{nn} = xiOS - zbEvalOnQuad1D{nn};
      elseif strcmp(interface, 'none')
        HOS{nn} = xiOSAlg(quadPhysPts1D{nn,1}, quadPhysPts1D{nn,2}, tOld) - zbEvalOnQuad1D{nn};
      else
        error('Invalid interface.');
      end % if
			
			assert(isequal(HOS{nn}<= 0, zeros(K*R1D,1)), 'HOS must be positive.');
			
    	uuH = cEdgeInt{2,nn}.^2 ./ HOS{nn};
      uvH = cEdgeInt{2,nn} .* cEdgeInt{3,nn} ./ HOS{nn};
      vvH = cEdgeInt{3,nn}.^2 ./ HOS{nn};
      quadraticH = 0.5 * gConst * HOS{nn}.^2;
      nonLinearity = nonLinearity + 0.5 * [ globROS{nn,1} * (uuH + quadraticH) + globROS{nn,2} * uvH; 
                                            globROS{nn,1} * uvH + globROS{nn,2} * (vvH +quadraticH) ];
			if strcmp(fluxType, 'Lax-Friedrichs')
				lambda = computeLaxFriedrichsCoefficientSWE('openSea', averaging, nn, [], kronNuE0T, [], [], [], [], [], HL, HOS, cEdgeInt, [], gConst);
				riemOS = riemOS + globVOS{nn} * (lambda .* (cEdgeInt{1,nn} - HOS{nn}));
			else
				error('Unknown type of flux approximation.');
			end % if
    elseif openSeaBdrs
      if strcmp(interface, 'ADCIRC')
        HOS{nn} = xiOS - zbEvalOnQuad1D{nn};
      elseif strcmp(interface, 'none')
        HOS{nn} = xiOSAlg(quadPhysPts1D{nn,1}, quadPhysPts1D{nn,2}, tOld) - zbEvalOnQuad1D{nn};
      else
        error('Invalid interface.');
      end % if
      uuH = cEdgeInt{2,nn}.^2 ./ HOS{nn};
      uvH = cEdgeInt{2,nn} .* cEdgeInt{3,nn} ./ HOS{nn};
      vvH = cEdgeInt{3,nn}.^2 ./ HOS{nn};
      quadraticH = 0.5 * gConst * HOS{nn}.^2;
			nonLinearity = nonLinearity + [ globROS{nn,1} * (uuH + quadraticH) + globROS{nn,2} * uvH; 
                                      globROS{nn,1} * uvH + globROS{nn,2} * (vvH + quadraticH) ];
		end % if
  end % for
	
	% building and solving the system
  sysV = [ globL{1} - globLRI{1}; globL{2} - globLRI{2}; globL{3} - globLRI{3} ];
	switch scheme
		case 'explicit'
			% calculate solution at current time via forard Euler scheme
      sysY = sysY + dt * ( sysW \ ( sysV - ( ( linTerms - [sparse(K*N,3*K*N); tidePot{1}, sparse(K*N,2*K*N); tidePot{2}, sparse(K*N,2*K*N)] ) ...
                                             * sysY + [ riemOS; nonLinearity + bottomFricPart ] + riem ) ) );
    case 'semi-implicit'
			% calculate solution at current time via semi-implicit Euler scheme
			sysY = (sysW + dt * ( linTerms - [sparse(K*N,3*K*N); tidePot{1}, sparse(K*N,2*K*N); tidePot{2}, sparse(K*N,2*K*N)] ) ) \ ...
             ( sysW * sysY + dt * ( sysV - ( [ riemOS; nonLinearity + bottomFricPart ] + riem ) ) );
		otherwise
			error('Invalid scheme.')  
	end % switch

  cDG(:,:,1)  = reshape(sysY(        1 :   K*N), N, K).';
  cDG(:,:,2)  = reshape(sysY(  K*N + 1 : 2*K*N), N, K).';
  cDG(:,:,3)  = reshape(sysY(2*K*N + 1 : 3*K*N), N, K).';
  % correction for unknown c1
	cDG(:,:,1) = applyMinValueExceedance2DataDisc(cDG(:,:,1), corrSys, nStep, minTol, 1000);
  
  % visualization
  if mod(nStep,output) == 0
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
  % save data in stations
  if useStations
		if ~(mod(nStep,output) == 0) % velocities are not available yet
			UDG(:,:,1) = projectQuotientDisc2DG(cDG(:,:,2), cDG(:,:,1), 2*p, refElemPhiPhi);
			UDG(:,:,2) = projectQuotientDisc2DG(cDG(:,:,3), cDG(:,:,1), 2*p, refElemPhiPhi);
		end % if
    uLagrange  = projectDataDisc2DataLagr(UDG(:,:,1)       );
    vLagrange  = projectDataDisc2DataLagr(UDG(:,:,2)       );
    xiLagrange = projectDataDisc2DataLagr(cDG(:,:,1) + zbDG);
		for stations = 1:NSTAE
			triangles = triE(stations);
			elevationInStations(stations, :) = sum(xiLagrange(triangles,:), 1) / length(triangles);
			elevationValues(nStep+1, stations) = extrapolateValue(g, triangles, coordE(stations,:), elevationInStations(stations, :));
		end % for
		for stations = 1:NSTAV
			triangles = triV(stations);
			velocityInStations(stations, :, 1) = sum(uLagrange(triangles,:), 1) / length(triangles);
			velocityInStations(stations, :, 2) = sum(vLagrange(triangles,:), 1) / length(triangles);
			velocityValues(nStep+1, stations, 1) = extrapolateValue(g, triangles, coordV(stations,:), velocityInStations(stations, :, 1));
			velocityValues(nStep+1, stations, 2) = extrapolateValue(g, triangles, coordV(stations,:), velocityInStations(stations, :, 2));
		end % for
  end % if 
  % waitbar
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
