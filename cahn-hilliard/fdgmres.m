function [x, error, total_iters] = fdgmres(f0, f, xc, params, P, xinit)
% GMRES linear equation solver for use in Newton-GMRES solver
%
% C. T. Kelley, July 24, 1994
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters] = fdgmres(f0, f, xc, params, xinit)
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%            params(3) (Optional) = reorthogonalization method
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%
%         xinit = initial iterate. xinit=0 is the default. This
%              is a reasonable choice unless restarted GMRES
%              will be used as the linear solver.
%
% Output: x=solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires givapp.m, dirder.m

%% Startwerte und Initialisierung
errtol = params(1); kmax = params(2); reorth = 1;
if length(params) == 3
  reorth = params(3);
end
% right side of linear equation for the step is -f0 if the
% default initial iterate is used
b = -f0;
n = length(b);
% Use zero vector as initial iterate for Newton step unless
% the calling routine has a better idea (useful for GMRES(m)).
x = zeros(n,1); r = b;
if nargin == 6
  x = xinit;
  %z = P(x);
  r = b - dirder(xc, x, f, f0, P);
end
h = zeros(kmax); v = zeros(n,kmax); c = zeros(kmax+1,1); s = zeros(kmax+1,1);
rho = norm(r); error = zeros(kmax,1);
g = rho * eye(kmax+1,1);
% test for termination on entry
total_iters = 0;
if(rho < errtol)
  disp(' early termination ')
  return
end
v(:,1) = r/rho; beta = rho; k = 0;
%% GMRES iteration
while((rho > errtol) && (k < kmax))
  k=k+1;
  %z = P(v(:,k));
  % call directional derivative function
  v(:,k+1) = dirder(xc, v(:,k), f, f0, P);
  normav = norm(v(:,k+1));
  
  % Modified Gram-Schmidt
  h(1:k,k) = v(:,1:k)'*v(:,k+1);
  v(:,k+1) = v(:,k+1) - v(:,1:k) * h(1 : k,k);
  h(k+1,k) = norm(v(:,k+1));
  normav2 = h(k+1,k);
  
  % reorthogonalize?
  if  (reorth == 1 && normav + .001*normav2 == normav) || reorth ==  3
    for j=1:k
      hr = v(:,j)' * v(:,k+1);
      h(j,k) = h(j,k)+hr;
      v(:,k+1) = v(:,k+1) - hr * v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
  end
  
  % watch out for happy breakdown
  if(h(k+1,k) ~= 0)
    v(:,k+1)=v(:,k + 1)/h(k + 1, k);
  else
    break
  end
  
  %
  %   Form and store the information for the new Givens rotation
  %
  if k > 1
    h(1:k,k)=givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
  end
  %
  %   Don't divide by zero if solution has  been found
  %
  nu=norm(h(k:k+1,k));
  if nu~=0
    c(k) = h(k,k)/nu;
    s(k) = -h(k + 1,k)/nu;
    h(k,k) = c(k) * h(k,k) - s(k) * h(k + 1,k);
    h(k+1,k) = 0;
    g(k:k+1) = givapp(c(k),s(k),g(k : k + 1),1);
  end
  % Update the residual norm
  rho = abs(g(k + 1));
  error(k) = rho;
end% end while
% At this point either k > kmax or rho < errtol.
% It's time to compute x and leave.
y = h(1:k,1:k) \ g(1:k);
error = error(1 : k);
total_iters = k;
x = x + v(:,1 : k)*y;
