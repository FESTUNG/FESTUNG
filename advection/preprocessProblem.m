function problemData = preprocessProblem(problemData)
%% Triangulation.
problemData.g = problemData.generateGridData(problemData.hmax);
if problemData.isVisGrid,  visualizeGrid(problemData.g);  end
%% Globally constant parameters.
problemData.K           = problemData.g.numT;  % number of triangles
problemData.N           = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
problemData.tau         = problemData.tEnd / problemData.numSteps;  % time step size

% [K x 3] arrays that mark local edges (E0T) or vertices (V0T) that are 
% interior or have a certain boundary type.
problemData.g.markE0Tint  = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);
problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...
                            problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD), :)); 

% Precompute some repeatedly evaluated fields                          
problemData.g = computeDerivedGridData(problemData.g);       
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
if problemData.isSlopeLim
  problemData.basesOnQuad = computeTaylorBasesV0T(problemData.g, problemData.N, problemData.basesOnQuad);
end % if
%% Computation of matrices on the reference triangle.
problemData.hatM              = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);
problemData.hatG              = integrateRefElemDphiPhiPhi(problemData.N, problemData.basesOnQuad);
problemData.hatRdiagOnQuad    = integrateRefEdgePhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
problemData.hatRoffdiagOnQuad = integrateRefEdgePhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
if problemData.isSlopeLim
  globMTaylor = assembleMatElemPhiTaylorPhiTaylor(problemData.g, problemData.N);
  problemData.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(problemData.g, problemData.N, problemData.basesOnQuad);
  problemData.globMCorr = spdiags(1./diag(globMTaylor), 0, problemData.K * problemData.N, problemData.K * problemData.N) * globMTaylor;
end % if
end % function