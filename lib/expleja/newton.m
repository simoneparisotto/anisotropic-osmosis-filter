function [y, err_est_norm, m] = newton (h, A, w, xi, dd, ...
                                        abs_tol, rel_tol, nnorm)
% Function file: [y, err_est_norm, m] = newton (h, A, w, xi, dd, ...
%                                       abs_tol, rel_tol, nnorm)
%
% Description: Perform Newton interpolation at real Leja points.
%
% Perform Newton interpolation. If dd are the divided differences of a
% function f, then y = p_m(h*A)*w, where m is the polynomial of degree
% m which interpolates f at the points xi(1:m+1). The maximum degree of
% interpolation is max_m, which corresponds to length (xi) - 1. The
% output degree m corresponds to the number of matrix-vector products.
%
% The error estimate works as follows:
%
% y_0 = d(1)*w
% y_k = y_{k-1}+d(k+1)*(h*A-xi(1)*I)*...*(h*A-xi(k)*I)*w, k = 1,...,m
%
% The estimate for y_k would be
%   err_est_norm=c1+c2
%               =||d(k+2)*(h*A-xi(1)*I)*(h*A-xi(k+1)*I)*w||+ ...
%                ||d(k+3)*(h*A-xi(1)*I)*...*(h*A-xi(k+2)*I)*w||,
% but since the arguments c1 and c2 are computed,
% y_{k+2} is computed as well. Therefore, the result is y_{m} with estimate
%   ||d(m)*(h*A-xi(1)*I)*...*(h*A-xi(m-1)*I)||+
%   ||d(m+1)*(h*A-xi(1)*I)*...*(h*A-xi(m)*I)||.
% The parameter nnorm selects the type of norm, must be valid parameter for
% NORM(*,nnorm).
% At least three points xi are required.
%
% This file is part of the expleja project.
% Authors: Marco Caliari and Peter Kandolf
% Date: September 23, 2016
  y = dd(1) * w;
  m = 1; % one matrix-vector product
  w = (A * w) * h - xi(1) * w;
  y = y + dd(m + 1) * w;
  c1 = norm (dd(m + 1) * w, nnorm);
  m = 2;
  w = (A * w) * h - xi(m) * w;
  y = y + dd(m + 1) * w;
  c2 = norm (dd(m + 1) * w, nnorm);
  err_est_norm = c1 + c2;
  if (err_est_norm <= max (rel_tol * norm (y, nnorm), abs_tol))
    return
  end
  for m = 3:length(dd) - 1
    w = (A * w) * h - xi(m) * w;
    y = y + dd(m + 1) * w;
    c1 = c2;
    c2 = norm (dd(m + 1) * w, nnorm);
    err_est_norm = c1 + c2;
    if (err_est_norm <= max (rel_tol * norm (y, nnorm), abs_tol))
      return
    end
  end
end
%!test % check output degree
%! m = 9;
%! A = hilb(20);
%! v = ones(20,1);
%! xi = (1:m+1);
%! dd = (-(m+1)/2:(m+1)/2-1);
%! [y, err_est_norm, m_test] = newton(1,A,v,xi,dd,1e-6,1e-6,inf);
%! assert(m_test, m)
%!test % check output degree
%! m = 9;
%! A = diag(-2:2);
%! v = ones(5,1);
%! xi = (-4:5);
%! dd = exptaylordd(xi);
%! [y, err_est_norm, m_test] = newton(1,A,v,xi,dd,1e-6,1e-6,inf);
%! assert(y, expm(A)*v, -eps)
%! assert(m_test, 8)
%!test % check interpolation
%! m = 20;
%! xi = 2 * cos((0:m) * pi / m); % Chebyshev points
%! dd = exptaylordd(xi);
%! A = rand(10)-1/2;
%! v = rand(10,1);
%! h = 0.1;
%! y = newton(h,A,v,xi,dd,1e-6,1e-6,inf);
%! assert(y,expm(h*A)*v,1e-6)
