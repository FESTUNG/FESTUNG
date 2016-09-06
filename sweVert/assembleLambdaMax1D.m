function [ lambdaMax1D ] = assembleLambdaMax1D( g , p , hDG , sysU , sysH , globIntV , globIntW , globR1D, gNum )

% lambdaMax1DR = 1.5 * abs(globIntV{1} * sysU) + 0.5 * ... 
%     sqrt( globIntW{1} * sysU + gNum * globR1D{1} * sysH );
% lambdaMax1DL = 1.5 * abs(globIntV{2} * sysU) + 0.5 * ... 
%     sqrt( globIntW{2} * sysU + gNum * globR1D{2} * sysH );
% 
% lambdaMax1D = 0.5 * (lambdaMax1DR(1:g.NX-1) + lambdaMax1DL(2:g.NX));

global gPsi1Dbnd

uLeft   = globIntV{1} * sysU;
uRight  = globIntV{2} * sysU;

uLeft   = uLeft(1:g.NX-1);
uRight  = uRight(2:g.NX);

hatU    = 0.5 * (uLeft + uRight);

hLeft   = zeros(g.NX, 1);
hRight  = zeros(g.NX, 1);

for i = 1 : p+1
    hLeft   = hLeft     + kron(hDG(:,i), gPsi1Dbnd(1,i)); % Left element
    hRight  = hRight    + kron(hDG(:,i), gPsi1Dbnd(2,i)); % Right element
end  % for

hLeft   = hLeft(1:g.NX-1);
hRight  = hRight(2:g.NX);

hatH = 0.5 * (hLeft + hRight);

hatH(hatH < 10^(-6)) = 10^(-6);

lambdaMax1D = ( 1.5 * abs(hatU) + 0.5 * sqrt(hatU.*hatU + 4 * gNum * hatH) ) ...
                ./ hatH;

end  % function