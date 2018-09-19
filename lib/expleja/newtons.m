function [y, err_est_norm, m] = newtons (h, A, w, xi, dd,...
                                       abs_tol, rel_tol, nnorm)
% Function file: [y, err_est_norm, m] = newtons (h, A, w, xi, dd, ...
%                                       abs_tol, rel_tol, nnorm)
%
% Description: Perform Newton interpolation at conjugate Leja points.
%
% If dd are the divided differences of a function f, real for real
% input, then y = p_m(h*A)*w, where m is the polynomial of degree m
% which interpolates f at the points xi(1:m+1).
% The maximum degree of interpolation is max_m, which corresponds to
% length (xi) - 1. The output degree m corresponds to the number of
% matrix-vector products.
%
% The error estimate works as follows:
% y_0 = d(1)*w
% y_2 = y_0+real(d(2))*h*A*w+d(3)*(h*A)*(h*A)*w
% y_k = y_{k-2}+real(d(k))*((h*A)^(k-1)+imag(xi(k-2))^2*(h*A)^(k-3))*w+...
%       d(k+1)*((h*A)^k+imag(xi(k-2))^2*(h*A)^(k-2)), k=4,...,m
%
% The estimate for y_k would be
%   err_est_norm = c1+c2
%          = ||real(d(k+2))*((h*A)^(k+1)+imag(xi(k))^2*(h*A)^(k-1))*w||+...
%            ||d(k+3)*((h*A)^(k+2)+imag(xi(k))^2*(h*A)^k)||,
% but since the arguments c1 and c2 are computed,
% y_{k+2} is computed as well. Therefore, the result is y_{m} with estimate
%   ||real(d(m))*((h*A)^(m-1)+...
%   imag(xi(m-2))^2*(h*A)^(m-3))*w||+||d(m+1)*((h*A)^m+...
%   imag(xi(m-2))^2*(h*A)^(m-2))||.
% The parameter nnorm selects the type of norm, must be valid parameter for
% NORM(*,nnorm).
% At least three points xi are required.
%
% This file is part of the expleja project.
% Authors: Marco Caliari and Peter Kandolf
% Date: September 23, 2016

%  dd(1:2:end) = real (dd(1:2:end)); % odd dd should be already real
  y = w * dd(1);
  m = 2; % two matrix-vector products
  w = (A * w) * h;
  wtilde = (A * w) * h;
  y = y + real (dd(m)) * w + real (dd(m + 1)) * wtilde;
  c1 = norm (real (dd(m)) * w, nnorm);
  c2 = norm (dd(m + 1) * wtilde, nnorm);
  err_est_norm = c1 + c2;
  if (err_est_norm <= max (rel_tol * norm (y, nnorm), abs_tol))
     return
  end
  for m = 4:2:length(dd) - 1
    w = (A * wtilde) * h + imag (xi(m - 2)) ^ 2 * w;
    wtilde = (A * w) * h;
    y = y + real (dd(m)) * w + real (dd(m + 1)) * wtilde;
    c1 = norm (real (dd(m)) * w, nnorm);
    c2 = norm (dd(m + 1) * wtilde, nnorm);
    err_est_norm = c1 + c2;
    if (err_est_norm <= max (rel_tol * norm (y, nnorm), abs_tol))
      return
    end
  end
end
%!test
%! m = 10;
%! xi = 1i * [0.00000, 2.00000, -2.00000, 1.15470, -1.15470, ...
%!       1.67989, -1.67989, 0.55213, -0.55213, 1.88604, -1.88604];
%! dd = exptaylordd(xi);
%! A = rand(4)-1/2;
%! v = rand(4,1);
%! h = 0.1;
%! y = newtons(h,A,v,xi,dd,1e-4,1e-4,inf);
%! assert(y,expm(h*A)*v,1e-4)
%!test % corresponds to the real part of complex Newton
%! m = 10;
%! xi = 1i * [0.00000, 2.00000, -2.00000, 1.15470, -1.15470, ...
%!       1.67989, -1.67989, 0.55213, -0.55213, 1.88604, -1.88604];
%! dd = exptaylordd(xi);
%! A = rand(4)-1/2;
%! v = rand(4,1);
%! h = 0.1;
%! ys = newtons(h,A,v,xi,dd,0,0,inf);
%! y = real (newton(h,A,v,xi,dd,0,0,inf));
%! assert(y,ys,eps)
%!test % corresponds to complex Newton
%! m = 10;
%! xi = 1i * [0.00000, 2.00000, -2.00000, 1.15470, -1.15470, ...
%!       1.67989, -1.67989, 0.55213, -0.55213, 1.88604, -1.88604];
%! dd = exptaylordd(xi);
%! A = rand(4)-1/2+1i*(rand(4)-1/2);
%! v = rand(4,1);
%! h = 0.1;
%! ys = newtons(h,A,v,xi,dd,0,0,inf);
%! y = newton(h,A,v,xi,dd,0,0,inf);
%! assert(y,ys,4*eps)
