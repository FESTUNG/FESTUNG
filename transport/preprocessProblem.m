function problemData = preprocessProblem(problemData)
%% Triangulation.
switch problemData.gridSource
  case 'square'
    problemData = setdefault(problemData, 'g', domainSquare(problemData.hmax));
    
    % Set edge types
    problemData.g.idE = zeros(problemData.g.numE,1);
    problemData.g.idE(problemData.g.baryE(:, 2) == 0) = 3; % south
    problemData.g.idE(problemData.g.baryE(:, 1) == 1) = 3; % east
    problemData.g.idE(problemData.g.baryE(:, 2) == 1) = 3; % north
    problemData.g.idE(problemData.g.baryE(:, 1) == 0) = 3; % west
    problemData.g.idE0T = problemData.g.idE(problemData.g.E0T);
  case 'hierarchical'
    problemData = setdefault(problemData, 'g', execin('swe/domainHierarchy', [0 100 100 0], [0 0 100 100], problemData.hmax, problemData.refinement));

    % Set edge types
    problemData.g.idE = zeros(problemData.g.numE,1);
    problemData.g.idE(problemData.g.baryE(:, 2) == 0) = 3; % south
    problemData.g.idE(problemData.g.baryE(:, 1) == 100) = 3; % east
    problemData.g.idE(problemData.g.baryE(:, 2) == 100) = 3; % north
    problemData.g.idE(problemData.g.baryE(:, 1) == 0) = 3; % west
    problemData.g.idE0T = problemData.g.idE(problemData.g.E0T);
  otherwise
    error('Invalid gridSource given.')
end % switch

if problemData.isVisGrid,  visualizeGrid(problemData.g);  end
%% Globally constant parameters.
problemData = setdefault(problemData, 'K', problemData.g.numT);  % number of triangles
problemData = setdefault(problemData, 'tau', problemData.tEnd / problemData.numSteps); % time step size
problemData.N = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
problemData = setdefault(problemData, 'velN', problemData.N); % number of local DOFs for velocity

problemData.g.markE0Tint  = problemData.g.idE0T == 0;        % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrN = zeros(problemData.g.numT,3);     % [K x 3] mark local edges on the Neumann boundary
problemData.g.markE0TbdrD = ~(problemData.g.markE0Tint | problemData.g.markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary
problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...  % [K x 3] mark local vertices on the Dirichlet boundary
                            problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD),:));
problemData.g = execin('transport/computeDerivedGridData', problemData.g);       % Precompute some repeatedly evaluated fields
%% Extract often used parameters
K = problemData.K;  % number of triangles
N = problemData.N;  % number of degrees of freedom
velN = problemData.velN; % number of degrees of freedom for velocity
p = problemData.p;  % Approximation order
velp = (sqrt(8*max(velN)+1)-3)/2;
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(max(N, velN), struct, sort(unique([2*p, 2*p+1, 2*velp, 2*velp+1])));
if any(cell2mat(problemData.isSlopeLim))
  problemData.basesOnQuad = computeTaylorBasesV0T(problemData.g, N, problemData.basesOnQuad);
end % if
%% Computation of matrices on the reference triangle.
problemData.hatM              = integrateRefElemPhiPhi(N, problemData.basesOnQuad);
problemData.hatG              = integrateRefElemDphiPhiPhi([N, N, problemData.velN], problemData.basesOnQuad);
problemData.hatRdiagOnQuad    = integrateRefEdgePhiIntPhiIntPerQuad(N, problemData.basesOnQuad);
problemData.hatRoffdiagOnQuad = integrateRefEdgePhiIntPhiExtPerQuad(N, problemData.basesOnQuad);
refElemPhiPerQuad = execin('swe/integrateRefElemPhiPerQuad', N, problemData.basesOnQuad);
%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globT = assembleMatElemPhiPhi(problemData.g, refElemPhiPerQuad);
if any(cell2mat(problemData.isSlopeLim))
  globMTaylor = assembleMatElemPhiTaylorPhiTaylor(problemData.g, N);
  problemData.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(problemData.g, N, problemData.basesOnQuad);
  problemData.globMCorr = spdiags(1./diag(globMTaylor), 0, K*N, K*N) * globMTaylor;
end % if
%% Function handle
problemData.swe_projectDataQ0T2DataDisc = getFunctionHandle('swe/projectDataQ0T2DataDisc');
end % function