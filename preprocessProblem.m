function problemData = preprocessProblem(problemData)

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
problemData.g.g1D = problemData.generateGrid1D(problemData.numElem(1), problemData.g);

%% Additional mesh data
% [K x 4] mark local edges that are horizontal or vertical
problemData.g.markE0Th = [ true(problemData.g.numT, 2), false(problemData.g.numT, 2) ];
problemData.g.markE0Tv = ~problemData.g.markE0Th;

% [K x 4] mark local edges that are interior or have a certain boundary type
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g); 
problemData.g.markE0Tbdr = ~problemData.g.markE0Tint; 
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrBot = problemData.generateMarkE0TbdrBot(problemData.g);
problemData.g.markE0TbdrTop = problemData.generateMarkE0TbdrTop(problemData.g);
problemData.g.markE0TbdrLand = problemData.generateMarkE0TbdrLand(problemData.g);
problemData.g.markE0TbdrOS = problemData.generateMarkE0TbdrOS(problemData.g);
problemData.g.markE0TbdrRiv = problemData.generateMarkE0TbdrRiv(problemData.g);
problemData.g.markE0TbdrRad = problemData.generateMarkE0TbdrRad(problemData.g);

% [barK x 2] mark local vertices that are interior or have a certain boundary type
problemData.g.g1D.markV0Tint = problemData.generateMarkV0T1Dint(problemData.g.g1D);
problemData.g.g1D.markV0Tbdr = ~problemData.g.g1D.markV0Tint;
problemData.g.g1D.markV0TbdrLand = problemData.generateMarkV0T1DbdrLand(problemData.g.g1D);
problemData.g.g1D.markV0TbdrOS = problemData.generateMarkV0T1DbdrOS(problemData.g.g1D);
problemData.g.g1D.markV0TbdrRiv = problemData.generateMarkV0T1DbdrRiv(problemData.g.g1D);
problemData.g.g1D.markV0TbdrRad = problemData.generateMarkV0T1DbdrRad(problemData.g.g1D);

% prescribed horizontal velocity
problemData.g.markE0TbdrU = problemData.g.markE0TbdrLand | problemData.g.markE0TbdrBot;
% prescribed vertical velocity
problemData.g.markE0TbdrW = problemData.g.markE0TbdrBot;
% prescribed water height
problemData.g.markE0TbdrH = problemData.g.markE0TbdrOS | problemData.g.markE0TbdrRiv;
% prescribed diffusion
problemData.g.markE0TbdrQ = problemData.g.markE0TbdrTop | problemData.g.markE0TbdrOS | problemData.g.markE0TbdrRad;
% prescribed flow rate
problemData.g.markE0TbdrUH = problemData.g.markE0TbdrRiv;
problemData.g.g1D.markV0TbdrUH = problemData.g.g1D.markV0TbdrRiv;

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs on trapezoidals
problemData.barN = problemData.p + 1;  % number of local DOFs on intervals
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Configuration output.
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('Running testcase "%s".\n', problemData.testcase);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d x %d (%d) trapezoids.\n', ...
        problemData.p, problemData.N, problemData.numElem(1), problemData.numElem(2), problemData.g.numT);
fprintf('%d time steps from t = %g to %g.\n', problemData.numSteps, problemData.t0, problemData.tEnd);
fprintf('-------------------------------------------------------------------------------------------\n');

%% Lookup table for basis function.
problemData.basesOnQuad1D = computeBasesOnQuad1D(problemData.p, struct, [problemData.qOrd, problemData.qOrd+1]);
problemData.basesOnQuad2D = computeBasesOnQuadTensorProduct(problemData.p, struct, [problemData.qOrd, problemData.qOrd+1]);

%% Computation of matrices on the reference element.
problemData.hatM = integrateRefElemTetraPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatG = integrateRefElemTetraDphiPhiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatH = integrateRefElemTetraDphiPhi(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatQdiag = integrateRefEdgeTetraPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatQoffdiag = integrateRefEdgeTetraPhiIntPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatRdiag = integrateRefEdgeTetraPhiIntPhiIntPhiInt(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatRoffdiag = integrateRefEdgeTetraPhiIntPhiExtPhiExt(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatSdiag = integrateRefEdgeTetraPhiIntPerQuad(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.hatQPerQuad = integrateRefEdgeTetraPhiIntPhiIntPerQuad(problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

problemData.barHatM = integrateRefElem1DPhiPhi(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barHatG = integrateRefElem1DDphiPhiPhiPerQuad(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barHatSdiag = integrateRefEdgeTetraPhi1DIntPerQuad(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barHatPdiag = computePhiIntPhiIntPhiIntV0T1D(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barHatPoffdiag = computePhiIntPhiExtPhiExtV0T1D(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

problemData.tildeHatH = integrateRefElemTetraDphiPhi1D([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.tildeHatQdiag = integrateRefEdgeTetraPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.tildeHatQoffdiag = integrateRefEdgeTetraPhiIntPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.tildeHatPdiag = integrateRefEdgeTetraPhiIntPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.tildeHatPoffdiag = integrateRefEdgeTetraPhiIntPhiExtPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);

%% Computation of time-independent 1D matrices.
problemData.barGlobM = assembleMatElemPhiPhi(problemData.g.g1D, problemData.barHatM);
end % function