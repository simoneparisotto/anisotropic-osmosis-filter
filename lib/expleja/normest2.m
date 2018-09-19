function [n, c] = normest2 (A, tol, iter_max, varargin)
%
% function n = normest2 (A)
% function n = normest2 (Afun)
% function n = normest2 (A, tol)
% function n = normest2 (A, tol, iter_max)
% function n = normest2 (Afun, tol, iter_max, P1, P2, ...)
% function [n, c] = normest2 (A, ...)
% function [n, c] = normest2 (A, ...)
%
% Description: Estimate the 2-norm of a matrix.
%
% Estimate the 2-norm of a matrix using the power series algorithm.
% n is the (under)estimate, c the number of iterations,
% that is the number of products of type (A' * A) * x.
% Afun is dessigned to work like in normest1.
% FEVAL(@AFUN,FLAG,X) for the following values of
%      FLAG       returns
%      'dim'      N
%      'real'     1 if A is real, 0 otherwise
%      'notransp' A*X
%      'transp'   A'*X
%
% Original references:
% N. Higham, Estimating the matrix $p$-norm, Numer. Math. 62, 539-555
% (1992)
% J. D. Dixon, Estimating extremal eigenvalues and condition numbers
% of matrices, SIAM J. Numer. Anal., 20(4), 812-814 (1983)
%
% This file is part of the expleja project.
% Authors: Marco Caliari and Peter Kandolf
% Date: September 23, 2016


% if A*x != 0, then y = A' * A * x != 0, otherways
% norm (A*x) = x' * A' * A * x = 0, contraddiction
% now, A*y != 0, otherways A' * A * y = 0 = (A' * A) ^ 2 * x would be
% zero. But A' * A is symmetric, therefore A' * A = P' * D * P, hence
% (A' * A) ^ 2 = P' * D ^ 2 * P. If (A' * A) ^ 2 * x =
% P' * D ^ 2 * P * x = 0, then D ^ 2 * P * x = 0, thus D * P * x = 0,
% hence P' * D * P * x = 0 = A' * A * x, contraddiction
% now, by induction the sequence is always different from zero. So, it
% is enough to check at the beginning only.
  if (isnumeric (A))
     Aisnum = true;
     size_of_A = size (A, 1);
  else
    Aisnum = false;
    size_of_A = A ('dim', [], varargin{:});
  end
  if ((nargin < 2) || (isempty (tol)))
     tol = 1e-6;
  end
  if (nargin < 3)
    iter_max = 100;
  end
  x = randn (size_of_A, 1);
  x = x / norm (x, 2);
  if (Aisnum)
    y = A * x;
  else
    y = A ('notransp', x, varargin{:});
  end
  n = norm (y, 2);
  if (n ~= 0)
    y = y / n;
    if (Aisnum)
      z = A' * y;
    else
      z = A ('transp', y, varargin{:});
    end
    c = 1;
    zn = norm (z, 2);
    if (zn <= z' * x) % see Higham, in particular, zn = 0
      return
    end
    x = z / zn;
  else
    %% n == 0 means
    %% 1) x in Ker(A) (small probability) or
    %% 2) A is null
    %% let's take another initial vector
    x = randn (size_of_A, 1);
    x = x / norm (x, 2);
    c = 0;
  end
  n_old = n;
  while (c < iter_max)
    if (Aisnum)
      y = A * x;
    else
      y = A ('notransp', x, varargin{:});
    end
    n = norm (y, 2);
    if (abs (n - n_old) <= tol * n)
      %% if A is null, n = n_old = 0 and the algorithm exits
      %% here.
      break
    end
    y = y / n;
    if (Aisnum)
      z = A' * y;
    else
      z = A ('transp', y, varargin{:});
    end
    c = c + 1;
    zn = norm (z, 2);
    if (zn <= z' * x)
       break
    end
    x = z / zn;
    n_old = n;
  end
end
%!test # here Matlab fails
%! ret = normest2(realmax);
%! assert(ret,realmax)
%!test # here Matlab fails
%! ret = normest2(realmin);
%! assert(ret,realmin)
%!test
%! ret = normest2(zeros(10));
%! assert(ret,0);
%!test
%! A = zeros(100);
%! A(1,1) = 1;
%! ret = normest2(A);
%! assert(ret,1);
%!test
%! A = 2*rand(10)-1;
%! ret = normest2(A);
%! assert(ret,norm(A,2),-1e-4);
%!test
%! A = 2*rand(10)-1+1i*(2*rand(10)-1);
%! ret = normest2(A);
%! assert(ret,norm(A,2),-1e-4);
%!test # here Matlab fails
%! A = toeplitz([-2,1,0,0]);
%! ret = normest2(A);
%! assert(ret,norm(A,2),-1e-4)
%!test # here Matlab fails
%! A = toeplitz([-2,1,0,0]);
%! [ret, c] = normest2(A, 1e-6, 3);
%! assert(c,3)
%!function z = afun_A_P (flag, x, A, m)
%! switch flag
%! case "dim"
%!   z = length (A);
%! case "real"
%!   z = isreal (A);
%! case "transp"
%!   z = x; for i = 1:m, z = A' * z;, endfor
%! case "notransp"
%!   z = x; for i = 1:m, z = A * z;, endfor
%! endswitch
%!endfunction
%!test # functional input
%! A = 2*rand(10)-1;
%! m = 2;
%! Afun = @(flag, x) afun_A_P (flag, x, A, m);
%! ret = normest2 (Afun);
%! assert (ret, norm (A ^ 2), 1e-4)
