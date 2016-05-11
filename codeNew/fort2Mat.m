function [ g, fcAlg, zbAlg, NOLIBF, NWP, bottomFric, NTIP, NRAMP, F, fT, ramping, xiRI, uRI, vRI, NBFR, xiOSX, xiOST, G, NDTVAR, DT, STATIM, RNDAY, ...
           ITRANS, CONVCR, H0, NSTAE, XEL, YEL, NSTAV, XEV, YEV, NHSTAR, NHSINC] = fort2Mat(name)
        
[ ICS, NOLIBF, NWP, NCOR, NTIP, NRAMP, G, NDTVAR, DT, STATIM, REFTIM, RNDAY, ITRANS, CONVCR, DRAMP, H0, SLAM0, SFEA0, TAU, CF, CORI, NTIF, TPK, AMIGT, ...
  ETRF, FFT, FACET, NBFR, AMIG, FF, FACE, NSTAE, XEL, YEL, NSTAV, XEV, YEV, NHSTAR, NHSINC] = fort15Read(['fort_' name '.15']);
% , NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE, NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV
[g, DP, AGRID]  = fort14Read(['fort_' name '.14'], ICS, SLAM0, SFEA0);
[g, NEDGES, NEDNO, NLEDN, NRAEDN, NRIEDN, ETRI, UNRI, UTRI, NSEDN, EMO, EFA] = fort17Read(['fort_' name '.17'], g, NBFR);
[g, NLEDN, NRAEDN, NRIEDN, NSEDN] = generateGridDataFromADCIRC(g, NEDGES, NEDNO, NLEDN, NRAEDN, NRIEDN, NSEDN);
g.idE = zeros(g.numE, 1);
g.idE(NLEDN ) = 1;
g.idE(NRAEDN) = 2;
g.idE(NRIEDN) = 3;
g.idE(NSEDN ) = 4;
g.idE0T = g.idE(g.E0T); % local edge IDs
fprintf(['Grid construction for ' AGRID 'successful.\n']);
% general parameters
SLAM = g.coordV(:,1)/6378206.4/cos(SFEA0) + SLAM0;
SFEA = g.coordV(:,2)/6378206.4;
RNDAY = RNDAY*86400;
% parameters in PDE
if NCOR
	fcAlg = @(x1,x2) evaluateFuncFromVertexValues(g,2*7.29212E-5*sin(SFEA),x1,x2);
else
	fcAlg = @(x1,x2) CORI*(x1==x1);
end % if
zbAlg = @(x1,x2) evaluateFuncFromVertexValues(g,-DP,x1,x2);
if NWP == 0
	if NOLIBF == 0
		bottomFric = TAU;
	elseif NOLIBF == 1
		bottomFric = CF;
	else
		error('Invalid type of bottom friction.');
	end % if
elseif NWP == 1
	error('There is no method to read bottom friction values varying in space.');
else
	error('Invalid type of bottom friction variation.');
end % if
% Newtonian tidal potential
fX = cell(2,1);
fT = cell(2,NTIF);
if NRAMP
	ramping = @(t_days) ((t_days-STATIM+REFTIM)/DRAMP)*(t_days < STATIM-REFTIM+DRAMP) + (t_days >= STATIM-REFTIM+DRAMP);
else
	ramping = @(t_days) 1;
end % if
F = [];
if NTIP
	F = cell(2,2,NTIF); % dim x trigonometric x component
	for n = 1:NTIF
		fT{1,n} = @(t)  cos(AMIGT(n)*t);
		fT{2,n} = @(t) -sin(AMIGT(n)*t);
		semi = round(0.00014/AMIGT(n));
		if semi == 1
			fX{1} = TPK(n) * FFT(n) * ETRF(n) * cos(SFEA).^2 .* cos(pi/180*(FACET(n) + 2*SLAM));
			fX{2} = TPK(n) * FFT(n) * ETRF(n) * cos(SFEA).^2 .* sin(pi/180*(FACET(n) + 2*SLAM));
		elseif semi == 2
			fX{1} = TPK(n) * FFT(n) * ETRF(n) * sin(2*SFEA)  .* cos(pi/180*(FACET(n) +   SLAM));
			fX{2} = TPK(n) * FFT(n) * ETRF(n) * sin(2*SFEA)  .* sin(pi/180*(FACET(n) +   SLAM));
		else
			error(['Unsupported frequency:' num2str(semi) '.']);
		end % if
		F{1,1,n} = G ./ (2*g.areaT) .* ( (fX{1}(g.V0T(:,2)) - fX{1}(g.V0T(:,1))).*g.B(:,2,2) - (fX{1}(g.V0T(:,3)) - fX{1}(g.V0T(:,1))).*g.B(:,2,1) );
		F{1,2,n} = G ./ (2*g.areaT) .* ( (fX{2}(g.V0T(:,2)) - fX{2}(g.V0T(:,1))).*g.B(:,2,2) - (fX{2}(g.V0T(:,3)) - fX{2}(g.V0T(:,1))).*g.B(:,2,1) );
		F{2,1,n} = G ./ (2*g.areaT) .* (-(fX{1}(g.V0T(:,2)) - fX{1}(g.V0T(:,1))).*g.B(:,1,2) + (fX{1}(g.V0T(:,3)) - fX{1}(g.V0T(:,1))).*g.B(:,1,1) );
		F{2,2,n} = G ./ (2*g.areaT) .* (-(fX{2}(g.V0T(:,2)) - fX{2}(g.V0T(:,1))).*g.B(:,1,2) + (fX{2}(g.V0T(:,3)) - fX{2}(g.V0T(:,1))).*g.B(:,1,1) );
	end
end % if
% River edges
if ~isequal(g.T0E(NRIEDN,2), zeros(length(NRIEDN),1))
	error('Some river edges are located in the interior of the domain.');
end % if
xiRI = zeros(g.numT, 1);
 uRI = zeros(g.numT, 1);
 vRI = zeros(g.numT, 1);
xiRI(g.T0E(NRIEDN,1)) = ETRI;
 uRI(g.T0E(NRIEDN,1)) = UNRI .* g.nuE(NRIEDN, 1) - UTRI .* g.nuE(NRIEDN, 2);
 vRI(g.T0E(NRIEDN,1)) = UNRI .* g.nuE(NRIEDN, 2) + UTRI .* g.nuE(NRIEDN, 1);
% Open sea edges
if ~isequal(g.T0E(NSEDN,2), zeros(length(NSEDN),1))
	error('Some open sea edges are located in the interior of the domain.');
end % if
xiOSX = cell(2,NBFR);
xiOST = cell(2,NBFR);
for i = 1:NBFR
	xiOST{1,i} = @(t)  cos(AMIG(i)*t);
	xiOST{2,i} = @(t) -sin(AMIG(i)*t);
	xiOSX{1,i} = sparse(g.numT,1);
	xiOSX{2,i} = sparse(g.numT,1);
	xiOSX{1,i}(g.T0E(NSEDN,1)) = FF(i) * EMO(i,:) .* cos(pi/180*(FACE(i) - EFA(i,:)));
	xiOSX{2,i}(g.T0E(NSEDN,1)) = FF(i) * EMO(i,:) .* sin(pi/180*(FACE(i) - EFA(i,:)));
end % if
end % function
