function problemData = preprocessProblem(problemData)

%% Function handles for time stepping
[~, ~, ~, subStepList] = getStepLists();
problemData.subStepHandles = getStepHandles(problemData.problemName, subStepList);

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs on trapezoidals
problemData.barN = problemData.p + 1;  % number of local DOFs on intervals
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size
problemData.isBottomFriction = ~strcmp(problemData.bottomFriction, 'none');

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
problemData.g.g1D = problemData.generateGrid1D(problemData.numElem(1), problemData.g);

%% Additional 2D mesh data: [K x 4] marker for local edges
% Edges that are horizontal or vertical
problemData.g.markE0Th = [ true(problemData.g.numT, 2), false(problemData.g.numT, 2) ];
problemData.g.markE0Tv = ~problemData.g.markE0Th;

% Interior edges
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g); 

% Vertical boundary edges
problemData.g.markE0Tbdr = sparse(~problemData.g.markE0Tint & problemData.g.markE0Tv); 

% Bottom boundary edges
problemData.g.markE0TbdrBot = sparse(problemData.generateMarkE0TbdrBot(problemData.g));

% Top boundary edges
problemData.g.markE0TbdrTop = sparse(problemData.generateMarkE0TbdrTop(problemData.g));

% Coupling boundary edges
problemData.g.markE0TbdrCoupling = sparse(problemData.generateMarkE0TbdrCoupling(problemData.g));

% Prescribed horizontal velocity
problemData.g.markE0TbdrU = sparse(problemData.generateMarkE0TbdrU(problemData.g));

% Prescribed water height
problemData.g.markE0TbdrH = sparse(problemData.generateMarkE0TbdrH(problemData.g));

% Prescribed diffusion
problemData.g.markE0TbdrQ = sparse(problemData.generateMarkE0TbdrQ(problemData.g));

% Boundary edges with Riemann solver applied
problemData.g.markE0TbdrRiem = sparse(problemData.generateMarkE0TbdrRiem(problemData.g));

assert(~any(problemData.g.markE0Th(:) & problemData.g.markE0TbdrRiem(:)), ...
       'Riemann solver for horizontal boundaries is not implemented!')
assert(nnz(problemData.g.markE0TbdrRiem & ~(problemData.g.markE0TbdrU | problemData.g.markE0TbdrH)) == 0, ...
       'Riemann solver specified for boundaries without Dirichlet data for U or H')

%% Additional 1D mesh data: [barK x 2] marker for local edges
% Interior vertices
problemData.g.g1D.markV0Tint = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0Tint(:, [4 3])) > 0;

% Boundary vertices
problemData.g.g1D.markV0Tbdr = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0Tbdr(:, [4 3])) > 0;

problemData.g.g1D.markV0TbdrU = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0TbdrU(:, [4 3])) > 0;
problemData.g.g1D.markV0TbdrH = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0TbdrH(:, [4 3])) > 0;

% Boundary vertices with Riemann solver applied
problemData.g.g1D.markV0TbdrRiem = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0TbdrRiem(:, [4 3])) > 0;

%% Function handles for problem-specific functions
problemData.fn_adaptFreeSurface = getFunctionHandle([problemData.problemName filesep 'adaptFreeSurface']);
problemData.fn_assembleMatEdgeTetraHorizPhiPhiNuBottomUp = getFunctionHandle([problemData.problemName filesep 'assembleMatEdgeTetraHorizPhiPhiNuBottomUp']);
problemData.fn_assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight']);
problemData.fn_assembleMatEdgeTetraPhiIntPhiIntFuncDisc1DIntNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatEdgeTetraPhiIntPhiIntFuncDisc1DIntNuHeight']);
problemData.fn_assembleMatElem1DDphiPhiFuncDiscHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatElem1DDphiPhiFuncDiscHeight']);
problemData.fn_assembleMatV0T1DPhiPhiFuncDiscNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatV0T1DPhiPhiFuncDiscNuHeight']);
problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight']);
problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu = getFunctionHandle([problemData.problemName filesep 'assembleVecEdgeTetraPhiIntFuncContHeightNu']);
problemData.fn_assembleVecEdgeTetraPhiIntFuncDiscIntHeightNu = getFunctionHandle([problemData.problemName filesep 'assembleVecEdgeTetraPhiIntFuncDiscIntHeightNu']);
problemData.fn_assembleVecV0T1DPhiIntFuncContNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleVecV0T1DPhiIntFuncContNuHeight']);
problemData.fn_assembleVecV0T1DPhiIntFuncDiscIntNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleVecV0T1DPhiIntFuncDiscIntNuHeight']);

%% Configuration output.
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('Running problem "%s" with testcase "%s".\n', problemData.problemName, problemData.testcase);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d x %d (%d) trapezoids.\n', ...
        problemData.p, problemData.N, problemData.numElem(1), problemData.numElem(2), problemData.g.numT);
fprintf('%d time steps from t = %g to %g.\n', problemData.numSteps, problemData.t0, problemData.tEnd);
fprintf('-------------------------------------------------------------------------------------------\n');

%% Lookup tables for basis function.
assert(problemData.qOrd <= problemData.qOrdMax, 'Maximum requested order of quadrature rule must be greater than or equal to qOrd')
problemData.basesOnQuad1D = computeBasesOnQuad1D(problemData.p, struct, problemData.qOrd : problemData.qOrdMax+1);
problemData.basesOnQuad2D = computeBasesOnQuadTensorProduct(problemData.p, struct, problemData.qOrd : problemData.qOrdMax+1);

%% Computation of matrices on the reference element.
problemData.hatM = integrateRefElemTetraPhiPhi(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatG = integrateRefElemTetraDphiPhiPhi(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatH = integrateRefElemTetraDphiPhi(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatQdiag = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatQoffdiag = integrateRefEdgePhiIntPhiExt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatRdiag = integrateRefEdgePhiIntPhiIntPhiInt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatSdiag = integrateRefEdgePhiIntPerQuad(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);

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

%% One-dimensional mass matrix in free-surface equation (XII)
problemData.barGlobM = assembleMatElemPhiPhi(problemData.g.g1D, problemData.barHatM);

%% Empty vectors and matrices for coupled problem
problemData.globJuCoupling = { sparse(problemData.g.numT * problemData.N, 1), sparse(problemData.g.numT * problemData.N, 1) };
problemData.globJwCoupling = sparse(problemData.g.numT * problemData.N, 1);
problemData.globJuuCoupling = { sparse(problemData.g.numT * problemData.N, 1), sparse(problemData.g.numT * problemData.N, 1) };
end % function