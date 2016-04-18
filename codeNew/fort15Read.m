function [ICS, NOLIBF, NWP, NCOR, NTIP, NRAMP, G, DT, STATIM, REFTIM, RNDAY, DRAMP, H0, SLAM0, SFEA0, TAU, CF, CORI, NTIF, TPK, AMIGT, ETRF, FFT, FACET, NBFR, AMIG, FF, FACE, NSTAE, XEL, YEL, NSTAV, XEV, YEV, NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE, NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV] = fort15Read(file)
fileID = fopen(file  );
RUNDES = fgets(fileID); % unused
RUNID  = fgets(fileID); % unused
IHOT   = cell2mat(textscan(fileID, '%f \n', 'CommentStyle', '!')); % not supported yet
OUTP   =				  textscan(fileID, '%s \n', 'CommentStyle', '!') ; % not supported yet
param  = cell2mat(textscan(fileID, '%f   ', 'CommentStyle', '!'));
ICS		 = param(1); assert(ICS		 == 1 || ICS		== 2, 'Invalid type of coordinate system.'				);
NOLIBF = param(2); assert(NOLIBF == 0 || NOLIBF == 1, 'Invalid type of bottom friction.'					); % test
NWP		 = param(3); assert(NWP		 == 0 || NWP    == 1,	'Invalid type of bottom friction variation.');
NCOR   = param(4); assert(NCOR   == 0 || NCOR   == 1, 'Invalid type of Coriolis parameter.'				);
NTIP   = param(5); assert(NTIP   == 0 || NTIP   == 1, 'Invalid type of Newtonian tide potential.' );
NWS    = param(6); assert(NWS    == 0 || NWS    == 1, 'Invalid type of wind stress.'							); % TODO
NRAMP  = param(7); assert(NRAMP  == 0 || NRAMP  == 1, 'Invalid type ramping.'											);
G			 = param(8); assert(isscalar(G),								'G has to be a real number.'								);
% This feature is not supported
NQUAD  = param(9);
XI		 = param(10:3:10+3*NQUAD-3);
YI		 = param(11:3:10+3*NQUAD-2);
W			 = param(12:3:10+3*NQUAD-1);
dataCountr = 10+3*NQUAD;
NDTVAR = param(dataCountr); dataCountr = dataCountr+1; assert(NDTVAR == 0 || NDTVAR == 1,		'Invalid type time increment selection.'								); % TODO
DT		 = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(DT) && DT > 0,				'Time increment has to be a positive real number.'			);
STATIM = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(STATIM),							'Start time has to be a real number.'										);
REFTIM = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(REFTIM),							'Reference time has to be a real number.'								); % TODO: usage?
RNDAY  = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(RNDAY) && RNDAY > 0, 'Invalid total length of simulation.'										);
IRK    = param(dataCountr); dataCountr = dataCountr+1; assert(IRK    == 0 || ...
																															IRK    == 1 || IRK    == 2,		'Invalid Runge-Kutta scheme.'														); % TODO 
ISLOPE = param(dataCountr); dataCountr = dataCountr+1; assert(ISLOPE == 0 || ISLOPE == 1 ...
																													 || ISLOPE == 2 || ISLOPE == 3,		'Invalid slope limiter or Riemann solver selection.'		); % TODO
ITRANS = param(dataCountr); dataCountr = dataCountr+1; assert(ITRANS == 0 || ITRANS == 1,		'The program must reach steady-state or run full time.' );
																					 assert(IRK >= 1 || ITRANS == 1,	'If no time-stepping scheme is used, the problem must be steady-state.' );
CONVCR = param(dataCountr); dataCountr = dataCountr+1; % TODO + assert
DRAMP  = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(DRAMP) && DRAMP >= 0, 'Ramping duration has to be positive.'									); % TODO: check UTBEST
H0		 = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(H0) && H0 > 0, 'Minimum cutoff depth has to be positive.'										);
SLAM0  = param(dataCountr); dataCountr = dataCountr+1;
SFEA0  = param(dataCountr); dataCountr = dataCountr+1; assert( isscalar(SLAM0) && isscalar(SFEA0), ['Reference coordinates for spherical ' ...
																																												'coordinates (CPP projection) must be real-valued numbers.']);
TAU		 = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(TAU), 'The linear bottom friction coefficient has to be a real number.'			);
CF		 = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(TAU), 'The non linear bottom friction coefficient has to be a real number.'	);
NVISC  = param(dataCountr); dataCountr = dataCountr+1; % TODO: check UTBEST
ESL		 = param(dataCountr); dataCountr = dataCountr+1; % TODO: check UTBEST
CORI	 = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(CORI), 'The Coriolis parameter has to be a real number.'											);
NTIF	 = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(NTIF) && round(NTIF) == NTIF, 'Invalid number of tide potential forcings.'		);
TIPOTAG = cell( NTIF,1);
TPK     = zeros(NTIF,1);
AMIGT   = zeros(NTIF,1);
ETRF    = zeros(NTIF,1);
FFT     = zeros(NTIF,1);
FACET   = zeros(NTIF,1);
if NTIF > 0
	for I = 1:NTIF
		TIPOTAG{I} = fgets(fileID);
		param			 = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
		TPK(I)		 = param(1);
		AMIGT(I)	 = param(2);
		ETRF(I)		 = param(3);
		FFT(I)		 = param(4);
		FACET(I)	 = param(5);
	end % for
	dataCountr = 6;
end % if
assert(isvector(TPK) && isvector(AMIGT) && isvector(ETRF) && isvector(FFT) && isvector(FACET),			 'Invalid tide potential parameters.');
NBFR   = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(NBFR) && round(NBFR) == NBFR, 'Invalid number of tidal forcings.' );
BOUNTAG = cell( NBFR,1);
AMIG	  = zeros(NBFR,1);
FF			= zeros(NBFR,1);
FACE		= zeros(NBFR,1);
if NBFR > 0
	for I = 1:NBFR
		BOUNTAG{I} = fgets(fileID);
		param			 = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
		AMIG(I)		 = param(1);
		FF(I)			 = param(2);
		FACE(I)		 = param(3);
	end % for
	dataCountr = 4;
end %if
assert(isvector(AMIG) && isvector(FF) && isvector(FACE), 'Invalid tidal forcing parameters.');
NOUTE   = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTSE  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTFE  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NSPOOLE = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NSTAE   = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(NSTAE) && round(NSTAE) == NSTAE, 'Invalid number of elevation recording stations');
if ICS == 1
	XEL  = param(dataCountr  :2:dataCountr+2*NSTAE-2);
	YEL  = param(dataCountr+1:2:dataCountr+2*NSTAE-1);
elseif ICS == 2
	SLEL = param(dataCountr  :2:dataCountr+2*NSTAE-2);
	SFEL = param(dataCountr+1:2:dataCountr+2*NSTAE-1);
	 XEL = 6378206.4*pi/180*(SLEL-SLAM0)*cos(pi/180*SFEA0);
	 YEL = 6378206.4*pi/180*SFEL;
else
  error('Invalid type of coordinate system.');
end % if
assert(isvector(XEL) && isvector(YEL), 'Invalid coordinates of some elevation recording stations.');
dataCountr = dataCountr+2*NSTAE;
NOUTV   = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTSV  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTFV  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NSPOOLV = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NSTAV   = param(dataCountr); dataCountr = dataCountr+1; assert(isscalar(NSTAV) && round(NSTAV) == NSTAV, 'Invalid number of velocity recording stations');
if ICS == 1
	XEV  = param(dataCountr  :2:dataCountr+2*NSTAV-2);
	YEV  = param(dataCountr+1:2:dataCountr+2*NSTAV-1);
else 
	SLEV = param(dataCountr	 :2:dataCountr+2*NSTAV-2);
	SFEV = param(dataCountr+1:2:dataCountr+2*NSTAV-1);
	 XEV = 6378206.4*pi/180*(SLEV-SLAM0)*cos(pi/180*SFEA0);
	 YEV = 6378206.4*pi/180*SFEV;
end % if
assert(isvector(XEV) && isvector(YEV), 'Invalid coordinates of some velocity recording stations.');
dataCountr = dataCountr+2*NSTAV;
NOUTGE   = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTSGE  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTFGE  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NSPOOLGE = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NOUTGV   = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTSGV  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
TOUTFGV  = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NSPOOLGV = param(dataCountr); dataCountr = dataCountr+1; % TODO, asserts
NHSTAR   = param(dataCountr); dataCountr = dataCountr+1; % TODO: check UTBEST
NHSINC   = param(dataCountr); % TODO: check UTBEST
fclose(fileID);
end % function
