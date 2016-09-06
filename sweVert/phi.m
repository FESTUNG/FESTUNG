function [ ret ] = phi( i , X , Y )

switch i
    case 1, ret = ones(size(X));
    case 2, ret = sqrt(3) * (1 - 2 * X);
    case 3, ret = sqrt(3) * (1 - 2 * Y);
    case 4, ret = 12 * X .* Y + 3 - 6 * ( X + Y );
    case 5, ret = sqrt(5) * ( (6 * X - 6) .* X + 1 );
    case 6, ret = sqrt(5) * ( (6 * Y - 6) .* Y + 1 );
    case 7, ret = sqrt(15) * ( 12 * X .* Y .* (X - 1) + 6 * X .* (1 - X) ...
                    + 2 * Y - 1 );
    case 8, ret = sqrt(15) * ( 12 * X .* Y .* (Y - 1) + 6 * Y .* (1 - Y) ...
                    + 2 * X - 1 );
    case 9, ret = 30 * ( 6 * X .* Y .* (X .* Y - X - Y + 1) + X .* X ...
                    + Y .* Y - X - Y + 1/6 );
end  % switch

end  % function