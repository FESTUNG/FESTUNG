function config = readConfigADCIRC(filename)
%% Open file for reading
assert(exist(filename, 'file') == 2, ['Config file "' filename '" does not exist!'])
fileID = fopen(filename, 'rt');

%% Read all parameters line-by-line and store them in a struct
config = struct;

% Text parameters
config.RUNDES = fgets(fileID);
config.RUNID = fgets(fileID);
config.IHOT = cell2mat(textscan(fileID, '%f \n', 'CommentStyle', '!'));
config.OUTP =	textscan(fileID, '%s \n', 'CommentStyle', '!');

% Numeric parameters
param  = cell2mat(textscan(fileID, '%f   ', 'CommentStyle', '!'));

% Coordinate system
config.ICS = param(1); 
assert(ismember(config.ICS, [1, 2]), 'Invalid type of coordinate system.');

% Non-linear bottom friction
config.NOLIBF = param(2); 
assert(ismember(config.NOLIBF, [0, 1]), 'Invalid type of bottom friction.');

% Spatially varying bottom friction
config.NWP = param(3); 
assert(ismember(config.NWP, [0, 1]), 'Invalid type of bottom friction variation.');

% Coriolis parameter
config.NCOR = param(4); 
assert(ismember(config.NCOR, [0, 1]), 'Invalid type of Coriolis parameter.');

% Newtonian tidal potential
config.NTIP = param(5); 
assert(ismember(config.NTIP, [0, 1]), 'Invalid type of Newtonian tide potential.');

% Wind stress
config.NWS = param(6);
assert(ismember(config.NWS, [0, 1]), 'Invalid type of wind stress.');

% Ramping parameters
config.NRAMP = param(7); 
assert(ismember(config.NRAMP, [0, 1]), 'Invalid type ramping.');

% Gravitational constant
config.G = param(8);
assert(isscalar(config.G), 'G has to be a real number.');

% Quadrature rules (Not used)
config.NQUAD = param(9);
config.XI = param(10:3:10+3*config.NQUAD-3);
config.YI = param(11:3:10+3*config.NQUAD-2);
config.W = param(12:3:10+3*config.NQUAD-1);
dataCountr = 10+3*config.NQUAD;

% Time step selection parameter
config.NDTVAR = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(ismember(config.NDTVAR, [0, 1]),'Invalid type time increment selection.');

% Time step size
config.DT = param(dataCountr);
dataCountr = dataCountr+1;
assert(isscalar(config.DT) && config.DT > 0, 'Time increment has to be a positive real number.');

% Start time in days
config.STATIM = param(dataCountr);
dataCountr = dataCountr+1;
assert(isscalar(config.STATIM), 'Start time has to be a real number.');

% Reference time in days
config.REFTIM = param(dataCountr);
dataCountr = dataCountr+1;
assert(isscalar(config.REFTIM), 'Reference time has to be a real number.');

% Simulation length in days
config.RNDAY = param(dataCountr);
dataCountr = dataCountr+1;
assert(isscalar(config.RNDAY) && config.RNDAY > 0, 'Invalid total length of simulation.');

% Runge-Kutta time stepping scheme
config.IRK = param(dataCountr);
dataCountr = dataCountr+1;
assert(ismember(config.IRK, [0, 1, 2]),	'Invalid Runge-Kutta scheme.');

% Slope reconstruction
config.ISLOPE = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(ismember(config.ISLOPE, [0, 1, 2, 3]), 'Invalid slope limiter or Riemann solver selection.' );

% Steady-state/transient simulation parameter
config.ITRANS = param(dataCountr); 
dataCountr = dataCountr+1;
assert(ismember(config.ITRANS, [0, 1]), 'The program must reach steady-state or run full time.');
assert(config.IRK >= 1 || config.ITRANS == 1,	'If no time-stepping scheme is used, the problem must be steady-state.');

% Convergence criteria 
config.CONVCR = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(config.ITRANS == 0 || config.CONVCR > 0, 'The tolerance for convergence of steady-state must be positive');

% Ramping period in days
config.DRAMP = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.DRAMP) && config.DRAMP >= 0, 'Ramping duration has to be positive.');

% Minimum allowed fluid-depth
config.H0 = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.H0) && config.H0 > 0, 'Minimum cutoff depth has to be positive.');

% Reference coordinates in degrees
config.SLAM0 = param(dataCountr); 
dataCountr = dataCountr+1;
config.SFEA0 = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.SLAM0) && isscalar(config.SFEA0), 'Reference coordinates for spherical coordinates (CPP projection) must be real-valued numbers.');

% Linear bottom friction coefficient
config.TAU = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.TAU), 'The linear bottom friction coefficient has to be a real number.');

% Non-linear bottom friction coefficient
config.CF = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.CF), 'The non linear bottom friction coefficient has to be a real number.');

% Lateral eddy viscosity (Not used)
config.NVISC = param(dataCountr); 
dataCountr = dataCountr+1; 
config.ESL = param(dataCountr); 
dataCountr = dataCountr+1; 

% Coriolis parameter
config.CORI	= param(dataCountr);
dataCountr = dataCountr+1;
assert(isscalar(config.CORI), 'The Coriolis parameter has to be a real number.');

% Tidal potential
config.NTIF = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.NTIF) && round(config.NTIF) == config.NTIF, 'Invalid number of tide potential forcings.');

config.TIPOTAG = cell(config.NTIF,1);
config.TPK = zeros(config.NTIF,1);
config.AMIGT = zeros(config.NTIF,1);
config.ETRF = zeros(config.NTIF,1);
config.FFT = zeros(config.NTIF,1);
config.FACET = zeros(config.NTIF,1);
if config.NTIF > 0
	for I = 1 : config.NTIF
		config.TIPOTAG{I} = fgets(fileID);
		param = cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
		config.TPK(I) = param(1);
		config.AMIGT(I) = param(2);
		config.ETRF(I) = param(3);
		config.FFT(I) = param(4);
		config.FACET(I) = param(5);
	end % for
	dataCountr = 6;
end % if
assert(isvector(config.TPK) && isvector(config.AMIGT) && isvector(config.ETRF) && isvector(config.FFT) && isvector(config.FACET), 'Invalid tide potential parameters.');

% Tidal forcing
config.NBFR = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.NBFR) && round(config.NBFR) == config.NBFR, 'Invalid number of tidal forcings.');
config.BOUNTAG = cell(config.NBFR,1);
config.AMIG = zeros(config.NBFR,1);
config.FF = zeros(config.NBFR,1);
config.FACE = zeros(config.NBFR,1);
if config.NBFR > 0
	for I = 1 : config.NBFR
		config.BOUNTAG{I} = fgets(fileID);
		param	= cell2mat(textscan(fileID, '%f', 'CommentStyle', '!'));
		config.AMIG(I) = param(1);
		config.FF(I) = param(2);
		config.FACE(I) = param(3);
	end % for
	dataCountr = 4;
end %if
assert(isvector(config.AMIG) && isvector(config.FF) && isvector(config.FACE), 'Invalid tidal forcing parameters.');

% Output elevation stations
config.NOUTE = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.TOUTSE = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.TOUTFE = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.NSPOOLE = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts

% Elevation recording stations
config.NSTAE = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.NSTAE) && round(config.NSTAE) == config.NSTAE, 'Invalid number of elevation recording stations');
switch config.ICS
  case 1
    config.XEL = param(dataCountr  :2:dataCountr+2*config.NSTAE-2);
    config.YEL = param(dataCountr+1:2:dataCountr+2*config.NSTAE-1);
  case 2
    % TODO put routine for CPP projection everywhere
    config.SLEL = param(dataCountr  :2:dataCountr+2*config.NSTAE-2);
    config.SFEL = param(dataCountr+1:2:dataCountr+2*config.NSTAE-1);
    config.XEL = 6378206.4*pi/180*(config.SLEL-config.SLAM0)*cos(pi/180*config.SFEA0);
    config.YEL = 6378206.4*pi/180*config.SFEL;
  otherwise
    error('Invalid type of coordinate system.');
end % switch
assert(isvector(config.XEL) && isvector(config.YEL), 'Invalid coordinates of some elevation recording stations.');
dataCountr = dataCountr+2*config.NSTAE;

% Output velocity stations
config.NOUTV = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.TOUTSV = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.TOUTFV = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.NSPOOLV = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts

% Velocity recording stations
config.NSTAV = param(dataCountr); 
dataCountr = dataCountr+1; 
assert(isscalar(config.NSTAV) && round(config.NSTAV) == config.NSTAV, 'Invalid number of velocity recording stations');
switch config.ICS
  case 1
    config.XEV = param(dataCountr  :2:dataCountr+2*config.NSTAV-2);
    config.YEV = param(dataCountr+1:2:dataCountr+2*config.NSTAV-1);
  case 2
    config.SLEV = param(dataCountr	 :2:dataCountr+2*config.NSTAV-2);
    config.SFEV = param(dataCountr+1:2:dataCountr+2*config.NSTAV-1);
    config.XEV = 6378206.4*pi/180*(config.SLEV-config.SLAM0)*cos(pi/180*config.SFEA0);
    config.YEV = 6378206.4*pi/180*config.SFEV;
  otherwise
    error('Invalid type of coordinate system.');
end % if
assert(isvector(config.XEV) && isvector(config.YEV), 'Invalid coordinates of some velocity recording stations.');
dataCountr = dataCountr+2*config.NSTAV;

% Output global elevation
config.NOUTGE = param(dataCountr); 
dataCountr = dataCountr+1; 
% TODO, asserts
config.TOUTSGE = param(dataCountr);
dataCountr = dataCountr+1; 
% TODO, asserts
config.TOUTFGE = param(dataCountr);
dataCountr = dataCountr+1;
% TODO, asserts
config.NSPOOLGE = param(dataCountr);
dataCountr = dataCountr+1;
% TODO, asserts

% Output global velocity
config.NOUTGV = param(dataCountr);
dataCountr = dataCountr+1;
% TODO, asserts
config.TOUTSGV = param(dataCountr);
dataCountr = dataCountr+1;
% TODO, asserts
config.TOUTFGV = param(dataCountr); 
dataCountr = dataCountr+1;
% TODO, asserts
config.NSPOOLGV = param(dataCountr);
dataCountr = dataCountr+1;
% TODO, asserts

% Output hot start
config.NHSTAR = param(dataCountr); 
dataCountr = dataCountr+1;
% TODO: check UTBEST
config.NHSINC = param(dataCountr);
% TODO: check UTBEST

%% Close file
fclose(fileID);
end % function
