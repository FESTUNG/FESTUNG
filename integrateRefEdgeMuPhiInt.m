%
function ret = integrateRefEdgeMuPhiInt(N, Nhybrid, basesOnQuad, basesOnGamma)
validateattributes(basesOnGamma, {'struct'}, {}, mfilename, 'basesOnGamma')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = Nhybrid-1;  qOrd = 2*p+1;  [~, W] = quadRule1D(qOrd);

ret = zeros(Nhybrid, N, 3, 2); % [N x N]
for nn = 1 : 3
    for i = 1 : Nhybrid
        for j = 1 : N
%             disp('====');
%             basesOnQuad.phi1D{qOrd}(:,i)
%             basesOnGamma.phi1D{qOrd}(:,j)
%             W'
%             W * basesOnQuad.phi1D{qOrd}(:,i, nn) .* basesOnGamma.phi1D{qOrd}(:,j)
%             sum( W' .* basesOnQuad.phi1D{qOrd}(:,i, nn) .* basesOnGamma.phi1D{qOrd}(:,j) )
            %       ret(i, j, nn) = sum( W' .* basesOnQuad.phi1D{qOrd}(:,i, nn) .* basesOnGamma.phi1D{qOrd}(:,j) );
            ret(i, j, nn, 1) =  W * (basesOnQuad.phi1D{qOrd}(:,j, nn) .* basesOnGamma.phi1D{qOrd}(:,i) );
            
            baseFlip = flipud(basesOnGamma.phi1D{qOrd});
            ret(i, j, nn, 2) =  W * (basesOnQuad.phi1D{qOrd}(:,j, nn) .* baseFlip(:,i) );
        end % for
    end % for
end % for
end % function
