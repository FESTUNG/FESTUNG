function problemData = preprocessProblem(problemData)
%% Triangulation.
if ~isfield(problemData, 'g'), problemData.g = domainSquare(problemData.hmax); end % if
% Alternative: ~isfield(problemData, 'g'), problemData.g = domainPolygon([0 1 1 0], [0 0 1 1], problemData.hmax); end % if
if problemData.isVisGrid,  visualizeGrid(problemData.g);  end
%% Globally constant parameters.
if ~isfield(problemData, 'K'), problemData.K           = problemData.g.numT; end % number of triangles
if ~isfield(problemData, 'tau'), problemData.tau         = problemData.tEnd / problemData.numSteps; end % time step size
problemData.N = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs

problemData.g.markE0Tint  = problemData.g.idE0T == 0;        % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrN = zeros(problemData.g.numT,3);     % [K x 3] mark local edges on the Neumann boundary
problemData.g.markE0TbdrD = ~(problemData.g.markE0Tint | problemData.g.markE0TbdrN); % [K x 3] mark local edges on the Dirichlet boundary
problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...  % [K x 3] mark local vertices on the Dirichlet boundary
                            problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD),:)); 
problemData.g = computeDerivedGridData(problemData.g);       % Precompute some repeatedly evaluated fields
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
if any(cell2mat(problemData.isSlopeLim))
  problemData.basesOnQuad = computeTaylorBasesV0T(problemData.g, problemData.N, problemData.basesOnQuad);
end % if
%% Computation of matrices on the reference triangle.
problemData.hatM              = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);
if isfield(problemData, 'velN')
  problemData.hatG              = integrateRefElemDphiPhiPhi([problemData.N, problemData.N, problemData.velN], problemData.basesOnQuad);
else
  problemData.hatG              = integrateRefElemDphiPhiPhi(problemData.N, problemData.basesOnQuad);
end % if
problemData.hatRdiagOnQuad    = integrateRefEdgePhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
problemData.hatRoffdiagOnQuad = integrateRefEdgePhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
if any(cell2mat(problemData.isSlopeLim))
  globMTaylor = assembleMatElemPhiTaylorPhiTaylor(problemData.g, problemData.N);
  problemData.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(problemData.g, problemData.N, problemData.basesOnQuad);
  problemData.globMCorr = spdiags(1./diag(globMTaylor), 0, problemData.K * problemData.N, problemData.K * problemData.N) * globMTaylor;
end % if
end % function