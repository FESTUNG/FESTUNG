function tableau = getDIRKtableau( order )

switch order
    case 1, tableau = getDIRK11Tableau;
    case 2, tableau = getDIRK22Tableau;
    case 3, tableau = getDIRK33Tableau;
    case 4, tableau = getDIRK54Tableau;
    otherwise
        msgID = 'getDIRKtableau:TimeIntegratorOrder';
        msg = 'Function is not defined for order>4.';
        baseException = MException(msgID,msg);
        causeID = 'getDIRKtableau:OrderTooLarge';
        causeMsg = sprintf('Value of order is %d', order);
        causeException = MException(causeID,causeMsg);
        baseException = addCause(baseException,causeException);
        throw(baseException);
end

end

function tableau = getDIRK11Tableau
lambda = 0.5 * (2. - sqrt(2.) );

A = [ 1 ];
B = A(end,:);
C = sum( A, 2);

%s
s = 1;
%Order of accuracy of the method
p = 1;
tableau = struct('A', A, 'B', B, 'C', C, 's', s, 'p', p);
end

function tableau = getDIRK22Tableau
lambda = 0.5 * (2. - sqrt(2.) );

A = [ lambda       ,  0;
      (1 - lambda) , lambda
    ];
B = A(end,:);
C = sum( A, 2);

%s
s = 2;
%Order of accuracy of the method
p = 2;
tableau = struct('A', A, 'B', B, 'C', C, 's', s, 'p', p);
end

function tableau = getDIRK33Tableau
alpha = 0.4358665215084589994160194511935568425292;
tau   = 0.5 * (1.+ alpha);
b1    = -0.25 * (6 * alpha * alpha - 16 * alpha + 1);
b2    =  0.25 * (6 * alpha * alpha - 20 * alpha + 5);

A = [ alpha       ,  0, 0;
      tau - alpha , alpha, 0
      b1, b2, alpha];
B = A(end,:);
C = sum( A, 2);

%s
s = 3;
%Order of accuracy of the method
p = 3;
tableau = struct('A', A, 'B', B, 'C', C, 's', s, 'p', p);
end

function tableau = getDIRK54Tableau
gamma = 1./4.; 
A = [     gamma,           0,        0,        0,     0;
          1./2.,       gamma,        0,        0,     0;
        17./50.,     -1./25.,    gamma,        0,     0;
     371./1360., -137./2720., 15./544.,    gamma,     0;
        25./24.,    -49./48., 125./16., -85./12., gamma
      ];
B = A(end,:);
C = sum( A, 2 );
%s
s = 5;
%Order of accuracy of the method
p = 4;
tableau = struct('A', A, 'B', B, 'C', C, 's', s, 'p', p);
end
