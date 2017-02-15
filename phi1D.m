function ret = phi1D(i, X)

%Mapping basis function from [-1,1] onto [0,1]
xi = 2 * ( X ) - 1;
% xi = 2 * ( X-1 )  1;

switch i
    case 1,  ret = ones( 1, size(X, 2) );
    case 2,  ret = xi .* sqrt(3);
    case 3,  ret = 0.5 * (3.*xi.^2 - 1) .* sqrt(5);
    case 4,  ret = 0.5 * ( 5 * xi.^3 - 3*xi ) .* sqrt(7);
    case 5,  ret = 0.125 * ( 35 * xi.^4 - 30*xi.^2 + 3) .* 3;
%     case 5,  ret = 0.125 * ( 63 * xi.^5 - 70*xi.^3 + 15*xi);
%     case 6,  ret = 1/16 * ( 231 * xi.^6 - 315*xi.^4 + 105*xi.^2 - 5);
    otherwise
        msgID = 'phi1D:BadPolynomialOrder';
        msg = 'Function is not defined for i>5 (p>4).';
        baseException = MException(msgID,msg);
        causeID = 'phi1D:iTooLarge';
        causeMsg = sprintf('Value of i is %d', i);
        causeException = MException(causeID,causeMsg);
        baseException = addCause(baseException,causeException);
        throw(baseException);
end % switch
end % function
