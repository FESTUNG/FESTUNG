function [errxi, erru, errv, erruH, errvH] = mainSWE(problem)
global gPhi2D gPhi1D gThetaPhi1D
refinement    = 0;
p             = 0;
OSRiem        = 1;
land          = 'Riemann';
scheme        = 'Runge-Kutta SSP/TVD'; % 'Runge-Kutta SSP/TVD' or 'semi-implicit Euler'
averaging     = 'full-harmonic';
wbarOpt       = 'off';
fileTypes     = 'vtk';
typeSlopeLim  = 'linear';       % Type of slope limiter (linear, hierarch_vert, strict)

fprintf('Read user input.\n');
[ interface, name, g, zbAlg, fcAlg, gConst, NDTVAR, NOLIBF, NWP, bottomFric, F0, F0Alg, rhsAlg, F1Alg, F2Alg, xiOSAlg, t0, tEnd, dt, numSteps, ...
	ordRK, ISLOPE, ITRANS, CONVCR, minTol, output, isVisParam, isVisu, isVisuH, isVisv, isVisvH, isVisxi, isSolAvail, xi, u, v, ...
  tidalDomain, F, fT, NRAMP, ramping, rhsRIAlg, xiRI, uRI, vRI, rhsOSAlg, NBFR, xiOSX, xiOST, useStations, NSTAE, triE, coordE, ...
  NSTAV, triV, coordV, NHSTAR, NHSINC ] = userInput(problem, refinement);

isSlopeLim = mod(ISLOPE, 2) == 1;
if ISLOPE == 0 || ISLOPE == 1
  fluxType = 'Roe'; % Not tested yet
else
  fluxType = 'Lax-Friedrichs';
end % if

%% Parameter check.
assert(ordRK >= 1 && ordRK <= 3, 'Order of Runge Kutta must be one to three.')
assert(~isSlopeLim || p > 0    , 'Slope limiting only available for p > 0.'   )

%% initial conditions
if isSolAvail % use exact solution at initial time
  xi0 = @(x1,x2) xi(x1,x2,t0);
   u0 = @(x1,x2)  u(x1,x2,t0);
   v0 = @(x1,x2)  v(x1,x2,t0);
else % cold start
  xi0 = @(x1,x2) zeros(size(x1));
  if problem == 4 % debug
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
markV0TbdrRI   = ismember(g.V0T, g.V0E(g.E0T(markE0TbdrRI),:)); % [K x 3] mark local vertices on the river boundary
markV0TbdrOS   = ismember(g.V0T, g.V0E(g.E0T(markE0TbdrOS),:)); % [K x 3] mark local vertices on the open sea boundary
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K);
g             = computeDerivedGridDataSWE(g, markE0Tint, markE0TbdrL, markE0TbdrRA, markE0TbdrRI, markE0TbdrOS);

landBdrs      = ~isequal(markE0TbdrL , zeros(K, 3));
radiationBdrs = ~isequal(markE0TbdrRA, zeros(K, 3));
riverBdrs     = ~isequal(markE0TbdrRI, zeros(K, 3));
openSeaBdrs   = ~isequal(markE0TbdrOS, zeros(K, 3));

[QX1, QX2, W] = quadRule2D(max(2*p,1)); R2D = length(W); % Determine the number of quadrature points on elements
zbEvalOnQuad2D = zbAlg(g.mapRef2Phy(1, QX1, QX2), g.mapRef2Phy(2, QX1, QX2));
zbEvalOnQuad2D = reshape(zbEvalOnQuad2D.', R2D*K, 1); 

[Q,W] = quadRule1D(max(2*p+1,1)); R1D = length(W); % Determine the number of quadrature points on edges

if strcmp(land, 'reflected')
  nuE0Tprod = cell(3,1);
  nuE0Tsqrd = cell(3,2);
  for nn = 1:3
    nuE0Tprod{nn}   = kron(g.nuE0T(:,nn,1) .* g.nuE0T(:,nn,2), ones(R1D,1));
    nuE0Tsqrd{nn,1} = kron(g.nuE0T(:,nn,1).^2, ones(R1D,1));
    nuE0Tsqrd{nn,2} = kron(g.nuE0T(:,nn,2).^2, ones(R1D,1));
  end % for
elseif strcmp(land, 'Riemann')
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
zbEvalOnQuadTheta1D = cell(3,3);                                           % This is used for boundary conditions
auxZb = cell(3,3);

if isVisParam
	zbLagrange = projectDataDisc2DataLagr(zbExact);
	visualizeDataLagr(g, zbLagrange, 'z_b', [name, '_z_b'], [], ['output_' name] , cd, fileTypes);
	fcLagrange = projectDataDisc2DataLagr(fcDG);
	visualizeDataLagr(g, fcLagrange, 'f_c', [name, '_f_c'], [], ['output_' name] , cd, fileTypes);
end % if

physQuadPts1D = cell(3,2);
kronNuE0T = cell(3,2);
for nn = 1:3
	[Q1, Q2] = gammaMap(nn, Q);
  aux = cell(2,1);
  for m = 1:2
    aux{m} = g.mapRef2Phy(m, Q1, Q2);
    physQuadPts1D{nn,m} = reshape(aux{m}.', R1D*K, 1);
    kronNuE0T{nn,m} = kron(g.nuE0T(:,nn,m), ones(R1D,1));
  end % for
	zbEvalOnQuad1D{nn} = zbAlg(aux{1}, aux{2});
  zbEvalOnQuad1D{nn} = reshape(zbEvalOnQuad1D{nn}.', R1D*K, 1);
  for np = 1:3
    [QP1, QP2] = theta(nn, np, Q1, Q2);
    zbEvalOnQuadTheta1D{nn,np} = zbAlg(g.mapRef2Phy(1, QP1, QP2), g.mapRef2Phy(2, QP1, QP2));
    zbEvalOnQuadTheta1D{nn,np} = reshape(zbEvalOnQuadTheta1D{nn,np}.', R1D*K, 1);
  end % for
end % for

if riverBdrs
%   [uRI, markE0TbdrRI, bsxfun(@times, markE0TbdrRI, uRI)]
%   [uRI markV0TbdrRI]
%   markV0TbdrRI
%   [uRI(markV0TbdrRI(:,1)) uRI(markV0TbdrRI(:,2)) uRI(markV0TbdrRI(:,3))]
%   [sum(uRI(g.V0T(:,1))) sum(uRI(g.V0T(:,2))) sum(uRI(g.V0T(:,3)))]
%   [sum(markV0TbdrRI(:,1)) sum(markV0TbdrRI(:,2)) sum(markV0TbdrRI(:,3))]
  if ~rhsRIAlg
    if isSlopeLim
      xiRIV0T = bsxfun(@times, markE0TbdrRI, xiRI); % TODO slopes
    end % if
    xiRI = kron(xiRI, ones(R1D,1));
     uRI = kron( uRI, ones(R1D,1));
     vRI = kron( vRI, ones(R1D,1));
  else
    if isSlopeLim
      xiRIV0T = computeFuncContV0T(g, @(x1, x2) xiRI(x1, x2, 0));
    end % if
  end % if  
end % if
if openSeaBdrs
  if ~rhsOSAlg
    if isSlopeLim
      xiOSV0T = cell(2, NBFR);
      xiOSV0TatStep = zeros(K, 3);
    end % for
    for i = 1 : NBFR
      if isSlopeLim
        for nn = 1:3
          xiOSV0T{1,i}(:,nn) = xiOSX{1,i} .* markE0TbdrOS(:,nn);
          xiOSV0T{2,i}(:,nn) = xiOSX{2,i} .* markE0TbdrOS(:,nn);
        end % for
        xiOSV0TatStep = xiOSV0TatStep + xiOST{1,i}(t0) * xiOSV0T{1,i} + xiOST{2,i}(t0) * xiOSV0T{2,i};
      end % if
      xiOSX{1,i} = kron(xiOSX{1,i}, ones(R1D,1));
      xiOSX{2,i} = kron(xiOSX{2,i}, ones(R1D,1));
    end % for
  else
    if isSlopeLim
      xiOSV0TatStep = computeFuncContV0T(g, @(x1, x2) xiOSAlg(x1, x2, 0));
    end % if
  end % if
end % if

% lookup table for basis functions
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

% well-balanced
refEdgePhiIntPhiExtPhiLin  = integrateRefEdgePhiIntPhiExtPhiLin (N);
refEdgePhiIntPhiIntPhiLin  = integrateRefEdgePhiIntPhiIntPhiLin (N);
refElemDphiPhiPhiLin       = integrateRefElemDphiPhiPhiLin      (N);
% refEdgePhiIntPhiLinPerQuad = integrateRefEdgePhiIntPhiLinPerQuad(N);

% only use zbDG for visualization and for error analysis
zbDG = projectFuncCont2DataDisc(g, @(x1,x2) zbAlg(x1,x2), 2*p, refElemPhiPhi);

for nn = 1:3 % TODO move this
  for np = 1:3
    auxZb{nn,np} = reshape(gThetaPhi1D{2*p+1}(:,:,nn,np) * zbDG.' * g.markE0TE0T{nn,np}.', R1D*K, 1);
  end % for
end % for

%% assembly of matrices corresponding to linear contributions
globQOS = assembleMatEdgePhiIntPhiIntNu    (g, markE0TbdrOS, refEdgePhiIntPhiInt,                      g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu    (g, markE0TbdrRA, refEdgePhiIntPhiInt,                      g.areaNuE0TbdrRA);
globH   = assembleMatElemDphiPhi           (g,               refElemDphiPhi                                            );
globM   = assembleMatElemPhiPhi						 (g,               refElemPhiPhi																						 );
globG   = assembleMatElemPhiPhiDfuncDiscLin(g,               refElemPhiPhiDphiLin, zbExact                             );
globD   = assembleMatElemPhiPhiFuncDiscLin (g,               refElemPhiPhiPhiLin,  fcDG                                );

% well-balanced
globOL  = assembleMatEdgePhiIntPhiIntPhiLinNu(g, markE0TbdrL , refEdgePhiIntPhiIntPhiLin, zbExact, g.areaNuE0TbdrL );
globORA = assembleMatEdgePhiIntPhiIntPhiLinNu(g, markE0TbdrRA, refEdgePhiIntPhiIntPhiLin, zbExact, g.areaNuE0TbdrRA);

% Roe - LF distinction
if strcmp(fluxType, 'Lax-Friedrichs')
  globQ = assembleMatEdgePhiPhiNu      (g, markE0Tint, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, g.areaNuE0Tint  );
  globO = assembleMatEdgePhiPhiPhiLinNu(g, markE0Tint, refEdgePhiIntPhiIntPhiLin, refEdgePhiIntPhiExtPhiLin, zbExact, g.areaNuE0Tint);
elseif strcmp(fluxType, 'Roe')
  jump = cell(3,1);
  globQ = assembleMatEdgePhiIntPhiIntNu      (g, markE0Tint, refEdgePhiIntPhiInt,                g.areaNuE0Tint);
  globO = assembleMatEdgePhiIntPhiIntPhiLinNu(g, markE0Tint, refEdgePhiIntPhiIntPhiLin, zbExact, g.areaNuE0Tint);
end % if

% well-balanced
globU   = assembleMatElemDphiPhiPhiLin(g, refElemDphiPhiPhiLin, zbExact);

sysW = blkdiag(globM, globM, globM);

if isSlopeLim
  computeTaylorBasesV0T(g, N);
  zbV0T = zeros(K,3);
  for nn = 1:3
    zbV0T(:,nn) = zbAlg(g.coordV0T(:, nn, 1), g.coordV0T(:, nn, 2));
  end % for
  globMTaylor     = assembleMatElemPhiTaylorPhiTaylor(g, N);
  globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(g, N);
  globMCorr       = spdiags(1./diag(globMTaylor), 0, K*N, K*N) * globMTaylor;
end % if

linTerms = [ sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1}, globQ{2} + globQOS{2} + globQRA{2} - globH{2}; 
             gConst * (globG{1} + globU{1} - globO{1} - globOL{1} - globORA{1}),	sparse(K*N,K*N), -globD;
             gConst * (globG{2} + globU{2} - globO{2} - globOL{2} - globORA{2}),	globD, sparse(K*N,K*N) ];

clear globD globG globH globQ globQOS globQRA globU globO globOL globORA

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
  if strcmp(land, 'natural')
    landFlux = sparse(K*N,1);
  elseif strcmp(land, 'Riemann')
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

globLRI    = cell(3,1);
globLRI{1} = sparse(K*N,1);
globLRI{2} = sparse(K*N,1);
globLRI{3} = sparse(K*N,1);

if riverBdrs && ~rhsRIAlg
  if ~NRAMP
    for nn = 1:3
        HRI = xiRI(:,nn) - zbEvalOnQuad1D{nn};
       uHRI =  uRI(:,nn) .*  HRI;
       vHRI =  vRI(:,nn) .*  HRI;
      uvHRI =  uRI(:,nn) .* vHRI;
      quadraticHRI = gConst * xiRI(:,nn) .* ( 0.5 * xiRI(:,nn) - zbEvalOnQuad1D{nn} );
      globLRI{1} = globLRI{1} + globRRI{nn,1} * uHRI + globRRI{nn,2} * vHRI;
      globLRI{2} = globLRI{2} + globRRI{nn,1} * (uRI(:,nn) .* uHRI + quadraticHRI) + globRRI{nn,2} * uvHRI;
      globLRI{3} = globLRI{3} + globRRI{nn,1} * uvHRI + globRRI{nn,2} * (vRI(:,nn) .* vHRI + quadraticHRI);
    end % for
%   else
%     HRI = sparse(K*N*R1D, 3);
%     uHRI = sparse(K*N*R1D, 3);
%     vHRI = sparse(K*N*R1D, 3);
  end % if
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

% system matrix for correction system 
corrSys = [ phi(1,0,0) phi(1,1,0) phi(1,0,1); ...
            phi(2,0,0) phi(2,1,0) phi(2,0,1); ...
            phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

% variable time step
if NDTVAR == 1
	dx = ( abs(g.coordV0T(:,1,1)-g.coordV0T(:,2,1)) + abs(g.coordV0T(:,2,1)-g.coordV0T(:,3,1)) + abs(g.coordV0T(:,3,1)-g.coordV0T(:,1,1)) ) / 3;
	dy = ( abs(g.coordV0T(:,1,2)-g.coordV0T(:,2,2)) + abs(g.coordV0T(:,2,2)-g.coordV0T(:,3,2)) + abs(g.coordV0T(:,3,2)-g.coordV0T(:,1,2)) ) / 3;
	avgDepth = -sum(zbAlg(g.coordV0T(:,:,1), g.coordV0T(:,:,2)), 2) / 3;
end % if

%% timestep zero
nStep =  0;
t     = t0;

% discrete data
cDG        = zeros(K,N,3);
cDG(:,:,1) = projectFuncCont2DataDisc(g, @(x1,x2)  xi0(x1,x2)                             , 2*p, refElemPhiPhi); % removed  - zbAlg(x1,x2)
cDG(:,:,2) = projectFuncCont2DataDisc(g, @(x1,x2) (xi0(x1,x2) - zbAlg(x1,x2)) .* u0(x1,x2), 2*p, refElemPhiPhi);
cDG(:,:,3) = projectFuncCont2DataDisc(g, @(x1,x2) (xi0(x1,x2) - zbAlg(x1,x2)) .* v0(x1,x2), 2*p, refElemPhiPhi);

if isSlopeLim
  if riverBdrs
    cDG(:,:,1) = applySlopeLimiterDisc(g, cDG(:,:,1), markV0TbdrRI, ramping(t0/86400) * xiRIV0T, globM, globMDiscTaylor, typeSlopeLim);
  end % if
  if openSeaBdrs
		cDG(:,:,1) = applySlopeLimiterDisc(g, cDG(:,:,1), markV0TbdrOS, xiOSV0TatStep, globM, globMDiscTaylor, typeSlopeLim);
  end % if
end % if

% correction for unknown c1
cDG(:,:,1) = applyMinValueExceedance2DataDisc(g, cDG(:,:,1) - zbDG, corrSys, nStep, minTol, 20); % put -zbDG here

cElemOnQuad = cell(3,1);
cEdgeIntOnQuad = cell(3,3);
cEdgeExtOnQuad = cell(3,3,3);
cEdgeExt2IntOnQuad = cell(3,3,3);

% visualize initial solution
UDG = zeros(K,N,2);
visNum = 0;
if isVisu || isVisv
  % compute velocities
  UDG = zeros(K, N, 2);
  UDG(:,:,1) = projectFuncDisc2DataDisc((cDG(:,:,2) * gPhi2D{max(2*p,1)}.') ./ (cDG(:,:,1) * gPhi2D{max(2*p,1)}.'), 2*p, refElemPhiPhi);
  UDG(:,:,2) = projectFuncDisc2DataDisc((cDG(:,:,3) * gPhi2D{max(2*p,1)}.') ./ (cDG(:,:,1) * gPhi2D{max(2*p,1)}.'), 2*p, refElemPhiPhi);
  if isVisu
    uLagrange = projectDataDisc2DataLagr(UDG(:,:,1)        );
    visualizeDataLagr(g,  uLagrange,  'u_h', [name '_u' ], visNum, ['output_' name] , cd, fileTypes);
  end % if
  if isVisuH
    uHLagrange = projectDataDisc2DataLagr(cDG(:,:,2)       );
    visualizeDataLagr(g, uHLagrange, 'uH_h', [name '_uH'], visNum, ['output_' name] , cd, fileTypes);
  end % if
  if isVisv
    vLagrange = projectDataDisc2DataLagr(UDG(:,:,2)        );
    visualizeDataLagr(g,  vLagrange,  'v_h', [name '_v' ], visNum, ['output_' name] , cd, fileTypes);
  end % if
  if isVisvH
    vHLagrange = projectDataDisc2DataLagr(cDG(:,:,3)       );
    visualizeDataLagr(g, vHLagrange, 'vH_h', [name '_vH'], visNum, ['output_' name] , cd, fileTypes);
  end % if
  if isVisxi
    xiLagrange = projectDataDisc2DataLagr(cDG(:,:,1) + zbDG);
    visualizeDataLagr(g, xiLagrange, 'xi_h', [name '_xi'], visNum, ['output_' name] , cd, fileTypes);
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
	if NDTVAR == 1
		dt = selectTimeStepSWE(dx, dy, avgDepth, gConst, dt, nStep);
	end % if
  nStep = nStep + 1;
  t     = t + dt;
  
  % TODO optimize this
  cDG(:,:,1) = cDG(:,:,1) + zbDG;
  
  % solution at old time step
  sysY = [ reshape(cDG(:,:,1).', K*N, 1); reshape(cDG(:,:,2).', K*N, 1); reshape(cDG(:,:,3).', K*N, 1) ];
	switch scheme
		case 'Runge-Kutta SSP/TVD'
			[timeLvls, omega] = rungeKuttaSSP(ordRK, dt, t - dt);
		  cDiscRK = cell(length(omega)+1, 1); cDiscRK{1} = sysY;
%       max(abs(sysY))
%       pause
		case 'semi-implicit Euler'
			timeLvls = t;
			cDiscRK = cell(2, 1); cDiscRK{1} = sysY;
		otherwise
			error('Invalid scheme.')
	end % switch
  
  %% Perform Runge-Kutta steps
  for rkStep = 1 : length(timeLvls)
  
		if NRAMP
			ramp = ramping(timeLvls(rkStep)/86400);
		end % if
    %% right hand sides / source terms / non solution dependent contributions
    if tidalDomain
      tidePot{1} = sparse(K*N,K*N);
      tidePot{2} = sparse(K*N,K*N);

      for n = 1 : size(F,3)
        tidePot{1} = tidePot{1} + fT{1,n}(timeLvls(rkStep))*F{1,1,n} + fT{2,n}(timeLvls(rkStep))*F{1,2,n};
        tidePot{2} = tidePot{2} + fT{1,n}(timeLvls(rkStep))*F{2,1,n} + fT{2,n}(timeLvls(rkStep))*F{2,2,n};
      end % for
      if NRAMP
        tidePot{1} = ramp * tidePot{1};
        tidePot{2} = ramp * tidePot{2};
      end % if
    end % if
    if rhsAlg
      if F0
        F0DG = projectFuncCont2DataDisc(g, @(x1,x2) F0Alg(x1,x2,timeLvls(rkStep)), 2*p, refElemPhiPhi);
        globL{1} = globM * reshape(F0DG', K*N, 1);
      end % if
      F1DG = projectFuncCont2DataDisc(g, @(x1,x2) F1Alg(x1,x2,timeLvls(rkStep)), 2*p, refElemPhiPhi);
      F2DG = projectFuncCont2DataDisc(g, @(x1,x2) F2Alg(x1,x2,timeLvls(rkStep)), 2*p, refElemPhiPhi);
      globL{2} = globM * reshape(F1DG', K*N, 1);
      globL{3} = globM * reshape(F2DG', K*N, 1);
    end % if
    
    if riverBdrs 
      if rhsRIAlg
        xiRIStep = xiRI(physQuadPts1D{nn,1}, physQuadPts1D{nn,2}, timeLvls(rkStep));
         vRIStep =  uRI(physQuadPts1D{nn,1}, physQuadPts1D{nn,2}, timeLvls(rkStep));
         uRIStep =  vRI(physQuadPts1D{nn,1}, physQuadPts1D{nn,2}, timeLvls(rkStep));
      elseif NRAMP
        xiRIStep = ramp * xiRI;
         uRIStep = ramp *  uRI;
         vRIStep = ramp *  vRI;
      end % if
    end % if
    
    if openSeaBdrs
      if strcmp(interface, 'ADCIRC') % TODO better cases
        xiOS = zeros(K*R1D, 1);
        for i = 1 : size(xiOSX,2)
          xiOS = xiOS + xiOST{1,i}(timeLvls(rkStep)) * xiOSX{1,i} + xiOST{2,i}(timeLvls(rkStep)) * xiOSX{2,i};
        end % for
        if NRAMP
          xiOS = ramp * xiOS;
        end % if
      else
        xiOS = xiOSAlg(physQuadPts1D{nn,1}, physQuadPts1D{nn,2}, timeLvls(rkStep));
      end % if
    end % if
    
    % initialize fields
    riem   = sparse(3*K*N,1);
    riemOS = sparse(K*N,1);
    if landBdrs && ~strcmp(land, 'natural')
      landFlux = sparse(K*N,1);
    end % if
    if riverBdrs && NRAMP % TODO keep in mind when using algebraic functions
      globLRI{1} = sparse(K*N,1);
      globLRI{2} = sparse(K*N,1);
      globLRI{3} = sparse(K*N,1);
    end % if
    
    %% solution dependent contributions
    % lookup table for solution
    cElemOnQuad{1} = reshape(gPhi2D{max(2*p,1)} * cDG(:,:,1).', R2D*K, 1);
    cElemOnQuad{2} = reshape(gPhi2D{max(2*p,1)} * cDG(:,:,2).', R2D*K, 1);
		cElemOnQuad{3} = reshape(gPhi2D{max(2*p,1)} * cDG(:,:,3).', R2D*K, 1);
    HElemOnQuad    = cElemOnQuad{1} - zbEvalOnQuad2D;
    
    %% element contributions
    % bottom friction part
    if NOLIBF == 0
      bottomFricPart = [globE * cDiscRK{rkStep}(K*N+1:2*K*N); globE * cDiscRK{rkStep}(2*K*N+1:end)];
    elseif NOLIBF == 1
      bottomFricNOLIBF = (cElemOnQuad{2}.^2 +cElemOnQuad{3}.^2).^0.5 ./ HElemOnQuad.^2;
      bottomFricPart = [globE * (bottomFricNOLIBF .* cElemOnQuad{2}); globE * (bottomFricNOLIBF .* cElemOnQuad{3})];
    else
      error('Invalid type of bottom friction.');
    end % if
    
    % nonlinearity
    uuH = cElemOnQuad{2}.^2 ./ HElemOnQuad;
    uvH = cElemOnQuad{2}.*cElemOnQuad{3} ./ HElemOnQuad;
    vvH = cElemOnQuad{3}.^2 ./ HElemOnQuad;
    quadraticH = 0.5 * gConst * cElemOnQuad{1}.^2;
    nonLinearity = [-globF{1} * (uuH + quadraticH)-globF{2} * uvH;
                    -globF{1} * uvH-globF{2} * (vvH + quadraticH)];
    
    % edge contributions
    for nn = 1 : 3
      cEdgeIntOnQuad{1,nn} = reshape(gPhi1D{2*p+1}(:,:,nn) * cDG(:,:,1).', R1D*K, 1);
			cEdgeIntOnQuad{2,nn} = reshape(gPhi1D{2*p+1}(:,:,nn) * cDG(:,:,2).', R1D*K, 1);
			cEdgeIntOnQuad{3,nn} = reshape(gPhi1D{2*p+1}(:,:,nn) * cDG(:,:,3).', R1D*K, 1);
      HEdgeIntOnQuad       = cEdgeIntOnQuad{1,nn} - zbEvalOnQuad1D{nn};
      if strcmp(averaging, 'full-harmonic') || strcmp(averaging, 'semi-harmonic')
        HL = sqrt(HEdgeIntOnQuad);
      end % if
      for np = 1:3
				aux = gThetaPhi1D{2*p+1}(:,:,nn,np) * cDG(:,:,1).';
				cEdgeExtOnQuad{1,nn,np} = reshape(aux, R1D*K, 1);
				cEdgeExt2IntOnQuad{1,nn,np} = reshape(aux * g.markE0TE0T{nn,np}.', R1D*K, 1);
				aux = gThetaPhi1D{2*p+1}(:,:,nn,np) * cDG(:,:,2).';
				cEdgeExtOnQuad{2,nn,np} = reshape(aux, R1D*K, 1);
				cEdgeExt2IntOnQuad{2,nn,np} = reshape(aux * g.markE0TE0T{nn,np}.', R1D*K, 1);
				aux = gThetaPhi1D{2*p+1}(:,:,nn,np) * cDG(:,:,3).';
				cEdgeExtOnQuad{3,nn,np} = reshape(aux, R1D*K, 1);
				cEdgeExt2IntOnQuad{3,nn,np} = reshape(aux * g.markE0TE0T{nn,np}.', R1D*K, 1);
        
        HEdgeExtOnQuad = cEdgeExtOnQuad{1,nn,np} - zbEvalOnQuadTheta1D{nn,np};
        HEdgeExt2IntOnQuad = cEdgeExt2IntOnQuad{1,nn,np} - auxZb{nn,np};
				
        uuH = cEdgeExtOnQuad{2,nn,np}.^2 ./ HEdgeExtOnQuad;
        uvH = cEdgeExtOnQuad{2,nn,np} .* cEdgeExtOnQuad{3,nn,np} ./ HEdgeExtOnQuad;
        vvH = cEdgeExtOnQuad{3,nn,np}.^2 ./ HEdgeExtOnQuad;
        quadraticH = 0.5*gConst*cEdgeExtOnQuad{1,nn,np}.^2;
        
        if strcmp(fluxType, 'Lax-Friedrichs') % TODO optimize this
          nonLinearity = nonLinearity + [ globRoffdiag{nn,np,1} * (uuH + quadraticH) + globRoffdiag{nn,np,2} * uvH;
                                          globRoffdiag{nn,np,1} * uvH + globRoffdiag{nn,np,2} * (vvH + quadraticH) ];
          lambda = computeLaxFriedrichsCoefficientSWE( 'interior', averaging, kronNuE0T, gConst, nn, cEdgeIntOnQuad, ...
                                                       HEdgeIntOnQuad, np, cEdgeExt2IntOnQuad, HEdgeExt2IntOnQuad, HL);
          riem = riem + [ globV{nn,np} * (lambda .* (cEdgeIntOnQuad{1,nn} - cEdgeExt2IntOnQuad{1,nn,np})); % changed for well-balanced
													globV{nn,np} * (lambda .* (cEdgeIntOnQuad{2,nn} - cEdgeExt2IntOnQuad{2,nn,np}));
													globV{nn,np} * (lambda .* (cEdgeIntOnQuad{3,nn} - cEdgeExt2IntOnQuad{3,nn,np})) ];
        elseif strcmp(fluxType, 'Roe')
          for I = 1:3
            jump{I} = cEdgeExt2IntOnQuad{I,nn,np} - cEdgeIntOnQuad{I,nn};
          end % for
          nonLinearity = nonLinearity + 2 * [ globRoffdiag{nn,1} * (uuH + quadraticH) + globRoffdiag{nn,2} * uvH;
                                              globRoffdiag{nn,1} * uvH + globRoffdiag{nn,2} * (vvH + quadraticH) ];
          [un, uHat, vHat, aHat] = computeRoe(); % TODO
          alphaLambda = min(un - aHat, 0) .* ( (aHat + un) .* jump{1} - kronNuE0T{nn,1} .* jump{2} - kronNuE0T{nn,2} .* jump{3} ) ./ (2 * a_hat);
          riem = riem + [ globV{nn,np} *   alphaLambda;
                          globV{nn,np} * ( alphaLambda .* (uHat - aHat .* kronNuE0T{nn,1}) );
                          globV{nn,np} * ( alphaLambda .* (vHat - aHat .* kronNuE0T{nn,2}) ) ];
          alphaLambda = min(un, 0) .* ( (uHat .* kronNuE0T{nn,2} - vHat .* kronNuE0T{nn,1}) .* jump{1} - kronNuE0T{nn,2} .* jump{2} + kronNuE0T{nn,1} .* jump{3} );
          riem = riem + [ sparse(K*R1D,1);
                          globV{nn,np} * (-alphaLambda .* kronNuE0T{nn,2});
                          globV{nn,np} * ( alphaLambda .* kronNuE0T{nn,1}) ];
          alphaLambda = min(un + aHat, 0) .* ( (aHat - un) .* jump{1} + kronNuE0T{nn,1} .* jump{2} + kronNuE0T{nn,2} .* jump{3} ) ./ (2 * a_hat);
          riem = riem + [ globV{nn,np} *   alphaLambda;
                          globV{nn,np} * ( alphaLambda .* (uHat + aHat .* kronNuE0T{nn,1}) );
                          globV{nn,np} * ( alphaLambda .* (vHat + aHat .* kronNuE0T{nn,2}) ) ];
        else
          error('Unknown type of flux approximation.');
        end % if
      end % for
      uuH = cEdgeIntOnQuad{2,nn}.^2 ./ HEdgeIntOnQuad;
      uvH = cEdgeIntOnQuad{2,nn} .* cEdgeIntOnQuad{3,nn} ./ HEdgeIntOnQuad;
      vvH = cEdgeIntOnQuad{3,nn}.^2 ./ HEdgeIntOnQuad;
      quadraticH = 0.5*gConst*cEdgeIntOnQuad{1,nn}.^2;
      nonLinearity = nonLinearity + [ globRdiag{nn,1} * (uuH + quadraticH) + globRdiag{nn,2} * uvH;
                                      globRdiag{nn,1} * uvH + globRdiag{nn,2} * (vvH + quadraticH) ];
      % land boundary contibutions
      if strcmp(land, 'natural')
        nonLinearity = nonLinearity + [globRL{nn,1}; globRL{nn,2}] * quadraticH;
      elseif strcmp(land, 'reflected')
        uHL = nuE0Tsqrd{nn,2} .* cEdgeIntOnQuad{2,nn} - nuE0Tprod{nn} .* cEdgeIntOnQuad{3,nn};
        vHL = nuE0Tsqrd{nn,1} .* cEdgeIntOnQuad{3,nn} - nuE0Tprod{nn} .* cEdgeIntOnQuad{2,nn};
        uuHL = uHL.^2 ./ HEdgeIntOnQuad;
        uvHL = uHL .* vHL ./ HEdgeIntOnQuad;
        vvHL = vHL.^2 ./ HEdgeIntOnQuad;
        landFlux = landFlux + globRL{nn,1} * uHL + globRL{nn,2} * vHL;
        nonLinearity = nonLinearity + [ globRL{nn,1} * (uuHL + quadraticH) + globRL{nn,2} * uvHL;
                                        globRL{nn,1} * uvHL + globRL{nn,2} * (vvHL + quadraticH) ];
      elseif strcmp(land, 'Riemann')
        uHR =  nuE0Tdiff{nn} .* cEdgeIntOnQuad{2,nn} - 2 * nuE0Tprod{nn} .* cEdgeIntOnQuad{3,nn};
        vHR = -nuE0Tdiff{nn} .* cEdgeIntOnQuad{3,nn} - 2 * nuE0Tprod{nn} .* cEdgeIntOnQuad{2,nn};
        if strcmp(fluxType, 'Lax-Friedrichs')
          lambda = computeLaxFriedrichsCoefficientSWE( 'land', [], kronNuE0T, gConst, nn, cEdgeIntOnQuad, [], [], [], [], HL, uHR, vHR );
          
          uuHR = uHR.^2 ./ HEdgeIntOnQuad;
          uvHR = uHR .* vHR ./ HEdgeIntOnQuad;
          vvHR = vHR.^2 ./ HEdgeIntOnQuad;
          
          landFlux = landFlux + 0.5 * (globRL{1} * (cEdgeIntOnQuad{2,nn} + uHR) + globRL{2} * (cEdgeIntOnQuad{3,nn} + vHR));
          nonLinearity = nonLinearity + 0.5 * ... % TODO schöner -> Roe
        [ globRL{nn,1} * (uuH + uuHR + 2 * quadraticH) + globRL{nn,2} * (uvH + uvHR) + globVL{nn} * (lambda .* (cEdgeIntOnQuad{2,nn} - uHR)); ...
          globRL{nn,1} * (uvH + uvHR) + globRL{nn,2} * (vvH + vvHR + 2 * quadraticH) + globVL{nn} * (lambda .* (cEdgeIntOnQuad{3,nn} - vHR)) ];
        else
          error('Unknown type of flux approximation.');
        end % if
      else
        error('Unknown type of land boundary discretization.')
      end % if
      
      % river boundary contributions
      if riverBdrs && NRAMP
          HRI = xiRIStep(:,nn) - zbEvalOnQuad1D{nn};
         uHRI =  uRIStep(:,nn) .*  HRI;
         vHRI =  vRIStep(:,nn) .*  HRI;
        uvHRI =  uRIStep(:,nn) .* vHRI;

        quadraticHRI = gConst * xiRIStep(:,nn) .* (0.5 * xiRIStep(:,nn) - zbEvalOnQuad1D{nn});
        globLRI{1} = globLRI{1} + globRRI{nn,1} * uHRI + globRRI{nn,2} * vHRI;
        globLRI{2} = globLRI{2} + globRRI{nn,1} * (uRIStep(:,nn) .* uHRI + quadraticHRI) + globRRI{nn,2} * uvHRI;
        globLRI{3} = globLRI{3} + globRRI{nn,1} * uvHRI + globRRI{nn,2} * (vRIStep(:,nn) .* vHRI + quadraticHRI);
        
        % well-balanced
%         globORI = assembleVecEdgePhiIntFuncDiscIntFuncContNu(g, markE0TbdrRI, refEdgePhiIntPhiLinPerQuad, zbExact, reshape(xiRIStep(:,nn), R1D, K).', g.areaNuE0TbdrRI);
%         nonLinearity = nonLinearity - gConst * [globORI{1}; globORI{2}];
        
      end % if
      
      % open sea boundary contributions
      if OSRiem
        quadraticH = quadraticH - gConst * cEdgeIntOnQuad{1,nn} .* zbEvalOnQuad1D{nn};
        nonLinearity = nonLinearity + 0.5 * [ globROS{nn,1} * (uuH + quadraticH) + globROS{nn,2} * uvH;
                                              globROS{nn,1} * uvH + globROS{nn,2} * (vvH + quadraticH) ];
        HOS{nn} = xiOS - zbEvalOnQuad1D{nn}; % TODO evtl HOS statt HOS{nn}
        assert(isequal(HOS{nn}<= 0, zeros(K*R1D,1)), 'HOS must be positive.');
        
        uuHOS = cEdgeIntOnQuad{2,nn}.^2 ./ HOS{nn};
        uvHOS = cEdgeIntOnQuad{2,nn} .* cEdgeIntOnQuad{3,nn} ./ HOS{nn};
        vvHOS = cEdgeIntOnQuad{3,nn}.^2 ./ HOS{nn};
        quadraticHOS = gConst * xiOS .* (0.5 * xiOS - zbEvalOnQuad1D{nn});
        nonLinearity = nonLinearity + 0.5 * [ globROS{nn,1} * (uuHOS + quadraticHOS) + globROS{nn,2} * uvHOS;
                                              globROS{nn,1} * uvHOS + globROS{nn,2} * (vvHOS + quadraticHOS) ];
        % well-balanced
%         globOOS = assembleVecEdgePhiIntFuncDiscIntFuncContNu(g, markE0TbdrOS, refEdgePhiIntPhiLinPerQuad, zbExact, reshape(xiOS, R1D, K).', g.areaNuE0TbdrOS);
%         nonLinearity = nonLinearity - gConst * [globOOS{1}; globOOS{2}]; % Note: This and same for river boundary are linear terms added here
        
        if strcmp(fluxType, 'Lax-Friedrichs')
          lambda = computeLaxFriedrichsCoefficientSWE( 'openSea', averaging, kronNuE0T, gConst, nn, cEdgeIntOnQuad, HEdgeIntOnQuad, ...
                                                       [], [], HEdgeExt2IntOnQuad, HL, [], [], HOS );
          riemOS = riemOS + globVOS{nn} * (lambda .* (HEdgeIntOnQuad - HOS{nn}));
        else
          error('Unknown type of flux approximation.');
        end % if
      elseif openSeaBdrs
        HOS{nn} = xiOS - zbEvalOnQuad1D{nn};
        assert(isequal(HOS{nn} <= 0, zeros(K*R1D,1)), 'HOS must be positive.');
        
        uuHOS = cEdgeIntOnQuad{2,nn}.^2 ./ HOS{nn};
        uvHOS = cEdgeIntOnQuad{2,nn} .* cEdgeIntOnQuad{3,nn} ./ HOS{nn};
        vvHOS = cEdgeIntOnQuad{3,nn}.^2 ./ HOS{nn};
        quadraticHOS = gConst * xiOS .* (0.5 * xiOS - zbEvalOnQuad1D{nn});
        nonLinearity = nonLinearity + [ globROS{nn,1} * (uuHOS + quadraticHOS) + globROS{nn,2} * uvHOS;
																				globROS{nn,1} * uvHOS + globROS{nn,2} * (vvHOS + quadraticHOS) ];
      end % if
    end % for
    
    % building and solving the system
    sysV = [ globL{1} - globLRI{1}; globL{2} - globLRI{2}; globL{3} - globLRI{3} ];
    switch scheme
      case 'Runge-Kutta SSP/TVD'
        % Computing the discrete time derivative
%         max(max(abs(linTerms)))
%         max(abs(nonLinearity))
%         max(abs(cDiscRK{rkStep}))
%         max(abs(riem))
%         max(abs(sysY))

        cDiscDot = sysW \ ( sysV - ( linTerms * cDiscRK{rkStep} - [sparse(K*N,K*N); tidePot{1}; tidePot{2}] * ...
                                      (sysY(1:K*N) - reshape(zbDG.', N*K,1)) ...
                           + [ riemOS; nonLinearity + bottomFricPart ] + riem ) ); % TODO tidePot
%                          max(abs(cDiscDot))
%                          pause
        if isSlopeLim
          for I = 1:1 % only limiting for height
            cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(K*N*(I-1)+1:K*N*I), [N K])', globM, globMDiscTaylor);
						
            cDiscDotTaylorLim = applySlopeLimiterTaylor(g, cDiscDotTaylor		, markV0TbdrRI, NaN(K,3), typeSlopeLim); % TODO hintereinander ok?
						cDiscDotTaylorLim = applySlopeLimiterTaylor(g, cDiscDotTaylorLim, markV0TbdrOS, NaN(K,3), typeSlopeLim); % TODO NaN's ok?
            cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
            cDiscDot(K*N*(I-1)+1:K*N*I) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', globM, globMDiscTaylor)', [K*N 1]);
          end % for
        end % if
        % Compute next step
%         max(abs(omega(rkStep) * sysY))
%         max(abs(dt * cDiscDot))
        
        cDiscRK{rkStep + 1} = omega(rkStep) * sysY + (1 - omega(rkStep)) * (cDiscRK{rkStep} + dt * cDiscDot);
        % Limiting the solution

% 				if isSlopeLim
%   if riverBdrs
%     if rhsRIAlg
%       xiRIV0T = computeFuncContV0T(g, @(x1, x2) xiRI(x1, x2, 0));
%     end % if
%     cDG(:,:,1) = applySlopeLimiterDisc(g, cDG(:,:,1), markV0TbdrRI, xiRIV0T, globM, globMDiscTaylor, typeSlopeLim) - zbDG;
%   end % if
%   if openSeaBdrs
%     if rhsOSAlg
%       xiSV0T = computeFuncContV0T(g, @(x1, x2) xiOSAlg(x1, x2, 0));
%     else
%       xiSV0T = zeros(K,3);
%       for i = 1 : size(xiOSX,2)
%         xiSV0T = xiSV0T + xiOST{1,i}(t0) * xiOSV0T{1,i} + xiOST{2,i}(t0) * xiOSV0T{2,i};
%       end % for
%     end % if
%     cDG(:,:,1) = applySlopeLimiterDisc(g, cDG(:,:,1), markV0TbdrOS, xiSV0T, globM, globMDiscTaylor, typeSlopeLim) - zbDG;
%   end % if
% end % if
        if isSlopeLim
          cDV0T = computeFuncContV0T(g, @(x1, x2) cDCont(t(rkStep), x1, x2));
          cDiscRK{rkStep + 1} = reshape(applySlopeLimiterDisc(g, reshape(cDiscRK{rkStep + 1}, [N K])', markV0TbdrD, cDV0T, globM, ...
                                            globMDiscTaylor, typeSlopeLim)', [K*N 1]);
        end % if
        
				if isSlopeLim
					if riverBdrs
						if rhsRIAlg
							xiRIV0T = computeFuncContV0T(g, @(x1, x2) xiRI(x1, x2, t(rkStep)));
            elseif NRAMP
							HRIV0T = xiRIV0T;
						end % if
						cDiscRK{rkStep + 1} = reshape(applySlopeLimiterDisc(g, reshape(cDiscRK{rkStep + 1}, [N K])', markV0TbdrRI, xiRIV0T, globM, ...
                                                                                  globMDiscTaylor, typeSlopeLim)', [K*N 1]);
					end % if
					if openSeaBdrs && false % TODO values of OSRiem
						if rhsOSAlg
							HOSV0T = computeFuncContV0T(g, @(x1, x2) xiOSAlg(x1, x2, t(rkStep))) - zbV0T;
						else
							HOSV0T = -zbV0T;
							for i = 1 : size(xiOSX,2)
								HOSV0T = HOSV0T + xiOST{1,i}(t0) * xiOSV0T{1,i} + xiOST{2,i}(t0) * xiOSV0T{2,i};
							end % for
						end % if
						cDiscRK{rkStep + 1} = reshape(applySlopeLimiterDisc(g, reshape(cDiscRK{rkStep + 1}, [N K])', markV0TbdrOS, HOSV0T, globM, ...
                                                                    globMDiscTaylor, typeSlopeLim)', [K*N 1]);
					end % if
				end % if
      case 'semi-implicit Euler' % no slope limiting, not well-balanced => will not work
        cDiscRK{rkStep + 1} = ( sysW + dt * ( linTerms - [sparse(K*N,3*K*N); tidePot{1}, sparse(K*N,2*K*N); ...
                                                                             tidePot{2}, sparse(K*N,2*K*N)] ) ) \ ...
																	( sysW * sysY + dt * ( sysV - ( [ riemOS; nonLinearity + bottomFricPart ] + riem ) ) );
      otherwise
        error('Invalid scheme.')
    end % switch
		cDG(:,:,1)  = reshape(cDiscRK{rkStep + 1}(        1 :   K*N), N, K).';
		cDG(:,:,2)  = reshape(cDiscRK{rkStep + 1}(  K*N + 1 : 2*K*N), N, K).';
		cDG(:,:,3)  = reshape(cDiscRK{rkStep + 1}(2*K*N + 1 : 3*K*N), N, K).';
		% correction for unknown c1
		cDG(:,:,1) = applyMinValueExceedance2DataDisc(g, cDG(:,:,1) - zbDG, corrSys, nStep, minTol, 20); % TODO error message in apply...
  end % for

	if ITRANS == 1
		l2 = sum((cDiscRK{end} - sysY).^2);
	end % if
  
  if mod(nStep, output) == 0 || (ITRANS == 1 && l2 < CONVCR)
%   if t >= 142560 && mod(nStep, 10) % debug gom
    UDG(:,:,1) = projectFuncDisc2DataDisc((cDG(:,:,2) * gPhi2D{max(2*p,1)}.') ./ (cDG(:,:,1) * gPhi2D{max(2*p,1)}.'), 2*p, refElemPhiPhi);
		UDG(:,:,2) = projectFuncDisc2DataDisc((cDG(:,:,3) * gPhi2D{max(2*p,1)}.') ./ (cDG(:,:,1) * gPhi2D{max(2*p,1)}.'), 2*p, refElemPhiPhi);
		visNum = visNum + 1;
		if isVisu
			uLagrange = projectDataDisc2DataLagr(UDG(:,:,1)        );
			visualizeDataLagr(g,  uLagrange,  'u_h', [name '_u' ] , visNum, ['output_' name] , cd, fileTypes);
		end % if
		if isVisuH
			uHLagrange = projectDataDisc2DataLagr(cDG(:,:,2)       );
			visualizeDataLagr(g, uHLagrange, 'uH_h', [name '_uH'], visNum, ['output_' name] , cd, fileTypes);
		end % if
		if isVisv
			vLagrange = projectDataDisc2DataLagr(UDG(:,:,2)        );
			visualizeDataLagr(g,  vLagrange,  'v_h', [name '_v' ] , visNum, ['output_' name] , cd, fileTypes);
		end % if
		if isVisvH
			vHLagrange = projectDataDisc2DataLagr(cDG(:,:,3)       );
			visualizeDataLagr(g, vHLagrange, 'vH_h', [name '_vH'], visNum, ['output_' name] , cd, fileTypes);
		end % if
		if isVisxi
			xiLagrange = projectDataDisc2DataLagr(cDG(:,:,1) + zbDG);
			visualizeDataLagr(g, xiLagrange, 'xi_h', [name '_xi'], visNum, ['output_' name] , cd, fileTypes);
		end % if
	end % if
  % save data in stations
  if useStations
		if ~(mod(nStep,output) == 0) % velocities are not available yet
			UDG(:,:,1) = projectFuncDisc2DataDisc((cDG(:,:,2) * gPhi2D{max(2*p,1)}.') ./ (cDG(:,:,1) * gPhi2D{max(2*p,1)}.'), 2*p, refElemPhiPhi);
			UDG(:,:,2) = projectFuncDisc2DataDisc((cDG(:,:,3) * gPhi2D{max(2*p,1)}.') ./ (cDG(:,:,1) * gPhi2D{max(2*p,1)}.'), 2*p, refElemPhiPhi);
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
	if NHSTAR == 1 && mod(nStep, NHSINC) == 0
		fprintf('\n\n')
		createHotStart([name  '_H_' ], cDG(:,:,1), t);
		createHotStart([name '_uH_' ], cDG(:,:,2), t);
		createHotStart([name '_vH_' ], cDG(:,:,3), t);
	end % if
  % waitbar
	if strcmp(wbarOpt, 'on')
	  percent = round( nStep / numSteps * 100);
		wbar = waitbar( percent/100, wbar, strcat( [ 'Time stepping:', ' ', num2str(percent), str ] ) );
	end % if
	if ITRANS == 1 && l2 < CONVCR
		fprintf('\n\nSteady state is reached.\n');
		break;
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
if useStations && NDTVAR == 0 % only supported for constant time step
	timeStations = t0:dt:tEnd;
	for stations = 1:NSTAE
		if stations ~= 1
			figure;
		end % if
		plot(timeStations, elevationValues(:, stations));
		title(['free surface elevation at station ', num2str(stations)]);
	end % for
	for stations = 1:NSTAV
		figure;
		plot(timeStations, velocityValues(:, stations, 1));
		title(['x-Component velocity at station ', num2str(stations)]);
		figure;
		plot(timeStations, velocityValues(:, stations, 2));
		title(['y-Component velocity at station ', num2str(stations)]);
	end % for
	cd(['output_' name])
	save('elevationInStations', 'elevationValues');
	save( 'velocityInStations',  'velocityValues');
	cd ..
end % if
end % function
