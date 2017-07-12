function [t, A, b, c] = rungeKuttaImplicit(ord, tau, t0)
switch ord
  case 1
    A = 1;
    
  case 2
    lambda = 0.5 * (2 - sqrt(2));
    A = [     lambda,      0 ; 
          (1-lambda), lambda ];
    
  case 3
    alpha = 0.4358665215084589994160194511935568425292;
    beta  = 0.5 * (1.+ alpha);
    b1    = -0.25 * (6 * alpha * alpha - 16 * alpha + 1);
    b2    =  0.25 * (6 * alpha * alpha - 20 * alpha + 5);

    A = [        alpha,     0,     0 ; 
          beta - alpha, alpha,     0 ; 
                    b1,    b2, alpha ];
    
  case 4
    gamma = 1./4.; 
    A = [     gamma,           0,        0,        0,     0;
              1./2.,       gamma,        0,        0,     0;
            17./50.,     -1./25.,    gamma,        0,     0;
         371./1360., -137./2720., 15./544.,    gamma,     0;
            25./24.,    -49./48., 125./16., -85./12., gamma ];
    
  otherwise
    msgID = 'rungeKuttaImplicit:TimeIntegratorOrder';
    msg = 'Function is not defined for order<1 or order>4.';
    baseException = MException(msgID,msg);
    causeID = 'rungeKuttaImplicit:OrderTooLarge';
    causeMsg = sprintf('Value of order is %d', ord);
    causeException = MException(causeID,causeMsg);
    baseException = addCause(baseException,causeException);
    throw(baseException);
end % switch
    
b = A(end, :);
c = sum(A, 2);
t = t0 + c * tau;
end % function
