function [ ret ] = gradPhi( i , m , X , Y )

switch m
    case 1
        switch i
            case 1, ret = zeros(size(X));
            case 2, ret = -2 * sqrt(3) * ones(size(X));
            case 3, ret = zeros(size(X));
            case 4, ret = 12 * Y - 6;
            case 5, ret = 6 * sqrt(5) * (2 * X - 1);
            case 6, ret = zeros(size(X));
            case 7, ret = 6 * sqrt(15) * (2 * X - 1) .* (2 * Y - 1);
            case 8, ret = 2 * sqrt(15) * (6 * Y .* Y - 6 * Y + 1);
            case 9, ret = 30 * (2 * X - 1) .* (6 * Y .* Y - 6 * Y + 1);
        end  % switch i
    case 2
        switch i
            case 1, ret = zeros(size(Y));
            case 2, ret = zeros(size(Y));
            case 3, ret = -2 * sqrt(3) * ones(size(Y));
            case 4, ret = 12 * X - 6;
            case 5, ret = zeros(size(Y));
            case 6, ret = 6 * sqrt(5) * (2 * Y - 1);
            case 7, ret = 2 * sqrt(15) * (6 * X .* X - 6 * X + 1);
            case 8, ret = 6 * sqrt(15) * (2 * X - 1) .* (2 * Y - 1);
            case 9, ret = 30 * (2 * Y - 1) .* (6 * X .* X - 6 * X + 1);
        end  % switch i
end  % switch m

end  % function

