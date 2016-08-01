function problemData = preprocessProblem(problemData)
%% Triangulation.
problemData = setdefault(problemData, 'g', domainSquare(problemData.hmax));
% problemData = setdefault(problemData, 'g',domainPolygon([0 1 1 0], [0 0 1 1], problemData.hmax)); % ALTERNATIVE
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
p = problemData.p;  % Approximation order
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(N, struct);
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
end % function