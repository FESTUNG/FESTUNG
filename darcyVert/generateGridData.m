function [ g ] = createDomain( width , NX , NZsub , NZsupra , intBound , upBound )

if size(intBound,2) ~= 1 || length(intBound) ~= NX+1
    error('Inner boundary is a vector of length NX+1.');
elseif size(intBound) ~= size(upBound)
    error('Boundaries must be of same size!');
end

%% Compute characteristic numbers
NZ       = NZsupra + NZsub;
g.numElem = [NX, NZ];
g.NX     = NX;
g.NZsupra = NZsupra;
g.numV   = (NX+1) * (NZ + 1);
g.numT   = NX * NZ;
g.numTsub = NX * NZsub;
g.numTsupra = NX * NZsupra;
g.numE   = NX * (NZ + 1) + NZ * (NX + 1);
g.deltaX = width / NX;
%% Create coordinate vector of vertices
g.coordV = zeros( g.numV , 2 );
g.coordV(:,1) = repmat( (0:g.deltaX:width)' , NZsub + NZsupra + 1 , 1 );
g.coordV(1:(NZsub+1)*(NX+1),2) = kron( (0 : NZsub)', intBound ./ NZsub );
g.coordV((NZsub+1)*(NX+1) + 1 :end,2) = kron(ones(NZsupra,1), intBound) ...
    + kron( (1 : NZsupra)', (upBound - intBound) ./ NZsupra );
%% Create Mapping: Trapezoid -> Vertices
g.V0T = zeros( g.numT , 4 );
g.V0T(:,1) = kron((0:NZsub+NZsupra-1)', (NX + 1) * ones(NX,1)) ...
    + repmat((1:NX)', NZsub + NZsupra , 1);
g.V0T(:,2) = g.V0T(:,1) + 1;
g.V0T(:,3) = g.V0T(:,2) + NX + 1;
g.V0T(:,4) = g.V0T(:,3) - 1;
%% Create Mapping: Edge -> Vertices
g.V0E = zeros( g.numE , 2 );
g.V0E(1:g.numT, :) = g.V0T(:,1:2);
g.V0E(g.numT:g.numT+NX,:) = g.V0T((NZ-1) * NX : end,4:-1:3);
g.V0E(NX * (NZ + 1) + 1 : g.numE - NZ,1) = g.V0T(:,1);
g.V0E(NX * (NZ + 1) + 1 : g.numE - NZ,2) = g.V0T(:,4);
g.V0E(g.numE - NZ + 1 : end , 1) = g.V0T(NX : NX : end, 2);
g.V0E(g.numE - NZ + 1 : end , 2) = g.V0T(NX : NX : end, 3);
%% Create Mapping: Trapezoid -> Edges
g.E0T = zeros( g.numT , 4 );
g.E0T(:,1) = ( 1 : g.numT );
g.E0T(:,2) = g.E0T(:,1) + NX;
g.E0T(:,3) = g.E0T(:,1) + NX * (NZ + 1) + 1;
g.E0T(:,4) = g.E0T(:,1) + NX * (NZ + 1);
g.E0T(NX : NX : end, 3) = (g.numE - NZ + 1 : g.numE);
%% Generate ID of edges
idE = zeros( g.numE, 1 );
idE(NZsub * NX + 1 : (NZsub + 1) * NX) = -1;
idE(1:NX) = 1;
idE(g.numE - NZ + 1 : end - NZsupra) = 2;
idE(g.numE - NZsupra + 1 : end) = 3;
idE(g.numT+1:g.numT+NX,:) = 4;
idE(g.numT + NX * (NZsub + 1) + 1 : NX :  g.numT + NX * (NZ + 1) , 1) = 5;
idE(g.numT + NX +1 : NX :  g.numT + NX * (NZsub + 1) , 1) = 6;
g.idE0T = idE(g.E0T);
g.idE = idE;
%% Create table of coordinates of vertices of Trapezodials
g.coordV0T = zeros( g.numT , 4 , 2 );
for i = 1 : 4
    g.coordV0T(:,i,:) = g.coordV( g.V0T(:,i) , : );
end  % for
g.coordV0Tsub = g.coordV0T(1:g.numTsub,:,:);
g.coordV0Tsupra = g.coordV0T(g.numTsub+1:g.numT,:,:);
%% Generate characteristic numbers for Jacobian Matrix
g.BA = zeros( g.numT , 1 );
g.DA = zeros( g.numT , 1 );
g.ACBD = zeros( g.numT , 1 );
g.BA = g.coordV0T(:,2,2) - g.coordV0T(:,1,2);
g.DA = g.coordV0T(:,4,2) - g.coordV0T(:,1,2);
g.CB = g.coordV0T(:,3,2) - g.coordV0T(:,2,2);
g.ACBD = g.coordV0T(:,1,2) + g.coordV0T(:,3,2) - g.coordV0T(:,2,2) ...
    - g.coordV0T(:,4,2);
g.BAsub = g.BA(1:g.numTsub);        g.BAsupra = g.BA(g.numTsub+1:g.numT);
g.DAsub = g.DA(1:g.numTsub);        g.DAsupra = g.DA(g.numTsub+1:g.numT);
g.ACBDsub = g.ACBD(1:g.numTsub);    g.ACBDsupra = g.ACBD(g.numTsub+1:g.numT);
                                    g.CBsupra = g.CB(g.numTsub+1:g.numT);
%% Generate Values of Nu * Length, where Nu is the outer unit normal
g.NuLengthE0T = zeros( g.numT , 4 , 2);
NuLengthE = ( g.coordV( g.V0E(:,2) , : ) - g.coordV( g.V0E(:,1) , : ) ) ...
    * [0, -1; 1, 0];
g.areaE = (g.coordV(g.V0E(:, 2), :) - g.coordV(g.V0E(:, 1), :));
g.nuE = NuLengthE ./ g.areaE;
g.NuLengthE0T(:,1,:) = NuLengthE( g.E0T(:, 1), : );
g.NuLengthE0T(:,2,:) = - NuLengthE( g.E0T(:, 2), : );
g.NuLengthE0T(:,3,:) = NuLengthE( g.E0T(:, 3), : );
g.NuLengthE0T(:,4,:) = - NuLengthE( g.E0T(:, 4), : );
g.NuLengthE0Tsub = g.NuLengthE0T(1:g.numTsub,:,:);
g.NuLengthE0Tsupra = g.NuLengthE0T(g.numTsub + 1 : g.numT,:,:);
%% Generate Length of Edge of Trapezodial
g.lengthE0T = zeros( g.numT , 4 );
g.lengthE0T(:,:) = sqrt( g.NuLengthE0T(:,:,1).^2 + g.NuLengthE0T(:,:,2).^2 );
g.lengthE0Tsub = g.lengthE0T(1:g.numTsub,:);
g.lengthE0Tsupra = g.lengthE0T(g.numTsub+1:g.numT,:);
%% Generate Mapping: {Edge} x Trapezodial -> Neighbouring Trapezodial
g.markE0TE0T = cell(4, 1);
g.markE0TE0Tsub = cell(4,1);
for i = 1 : 4
    g.markE0TE0T{i} = sparse( g.numT , g.numT );
end  % for
g.markE0TE0T{1}( sub2ind([g.numT,g.numT], NX+1:g.numT , 1:g.numT-NX) ) = 1;
g.markE0TE0T{2}( sub2ind([g.numT,g.numT], 1:g.numT-NX , NX+1:g.numT) ) = 1;
g.markE0TE0T{3}( sub2ind([g.numT,g.numT], 1:g.numT-1 , 2:g.numT) ) = 1;
g.markE0TE0T{3}( sub2ind([g.numT,g.numT], NX:NX:g.numT-NX , NX+1:NX:g.numT-NX+1) ) = 0;
g.markE0TE0T{4}( sub2ind([g.numT,g.numT], 2:g.numT , 1:g.numT-1) ) = 1;
g.markE0TE0T{4}( sub2ind([g.numT,g.numT], NX+1:NX:g.numT-NX+1 , NX:NX:g.numT-NX) ) = 0;
g.markE0TE0Tsub{1} = g.markE0TE0T{1}(1:g.numTsub,1:g.numTsub);
g.markE0TE0Tsub{2} = g.markE0TE0T{2}(1:g.numTsub,1:g.numTsub);
g.markE0TE0Tsub{3} = g.markE0TE0T{3}(1:g.numTsub,1:g.numTsub);
g.markE0TE0Tsub{4} = g.markE0TE0T{4}(1:g.numTsub,1:g.numTsub);
g.markE0TE0Tsupra{1} = g.markE0TE0T{1}(g.numTsub+1:g.numT,g.numTsub+1:g.numT);
g.markE0TE0Tsupra{2} = g.markE0TE0T{2}(g.numTsub+1:g.numT,g.numTsub+1:g.numT);
g.markE0TE0Tsupra{3} = g.markE0TE0T{3}(g.numTsub+1:g.numT,g.numTsub+1:g.numT);
g.markE0TE0Tsupra{4} = g.markE0TE0T{4}(g.numTsub+1:g.numT,g.numTsub+1:g.numT);
%% Generate Vector of areas of Trapezoids
g.areaT = zeros( g.numT , 1 );
g.areaT = 0.5 * g.deltaX * ( g.DA + g.coordV0T(:,3,2) - g.coordV0T(:,2,2) );
g.areaTsub = g.areaT(1:g.numTsub);
g.areaTsupra = g.areaT(g.numTsub + 1 : g.numT);
%% Generate Vector of Barycenters of Edges of Tapezoids
g.baryE0T = zeros( g.numT , 2 );
g.baryE = 0.5 * (g.coordV(g.V0E(:, 1), :) + g.coordV(g.V0E(:, 2), :));
for i = 1 : 4
    g.baryE0T(:, i, 1) = g.baryE(g.E0T(:,i), 1)';
    g.baryE0T(:, i, 2) = g.baryE(g.E0T(:,i), 2)';
end  % for
g.baryE0Tsub = g.baryE0T(1:g.numTsub,:,:);
g.baryE0Tsupra = g.baryE0T(g.numTsub+1:g.numT,:,:);
%% Generate Vector of Barycenters of Trapezoids
g.baryT = ( g.coordV0T(:,1,:) + g.coordV0T(:,2,:) + g.coordV0T(:,3,:) ...
    + g.coordV0T(:,4,:) ) / 4;
g.baryTsub = g.baryT(1:g.numTsub,:,:);
%% Generate Data for 1D Domain Interior Boundary
g.markE0Tint1D = ones( NX , 2 );
g.markE0TE0T1D = cell(2 , 1);
g.markE0TE0T1D{1} = sparse( NX , NX );
g.markE0TE0T1D{2} = sparse( NX , NX );
g.coordL1D = (0 : g.deltaX : NX*g.deltaX - 0.1*g.deltaX);
g.markE0Tint1D(1,2) = 0;
g.markE0Tint1D(NX,1) = 0;
g.markE0TE0T1D{1}( sub2ind([NX,NX], 1:NX-1, 2:NX) ) = 1;
g.markE0TE0T1D{2}( sub2ind([NX,NX], 2:NX, 1:NX-1) ) = 1;
end  % function