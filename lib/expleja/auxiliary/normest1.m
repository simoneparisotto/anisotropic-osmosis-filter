## Copyright (C) 2016 Marco Caliari
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{c} =} normest1 (@var{a})
## @deftypefnx {Function File} {@var{c} =} normest1 (@var{a}, @var{t})
## @deftypefnx {Function File} {@var{c} =} normest1 (@var{a}, @var{t}, @var{x0})
## @deftypefnx {Function File} {@var{c} =} normest1 (@var{afun}, @var{t}, @var{x0}, @var{p1}, @var{p2}, @dots{})
## @deftypefnx  {Function File} {[@var{c}, @var{v}] =} normest1 (@var{a}, @dots{})
## @deftypefnx  {Function File} {[@var{c}, @var{v}, @var{w}] =} normest1 (@var{a}, @dots{})
## @deftypefnx  {Function File} {[@var{c}, @var{v}, @var{w}, @var{it}] =} normest1 (@var{a}, @dots{})
##
## Estimate the 1-norm of the matrix A using a block algorithm.
##
## For a medium size matrix @var{A}, @code{norm (A, 1)} should be
## used instead. For a large sparse matrix, when only an estimate of the norm
## is needed, @code{normest1 (A)} might be faster. Moreover, it can be
## used for the estimate of the 1-norm of a linear
## operator @var{A} when matrix-vector products @code{A * X} and
## @code{A' * X} can be cheaply computed. In this case, instead of the
## matrix @var{A}, a function @code{AFUN (FLAG, X)} can be
## used. It should returns
##
## the dimension @var{N} of @var{A}, if @var{flag} is 'dim'
##
## true if @var{A} is a real operator, if @var{flag} is 'real'
##
## the result @code{A * X}, if @var{flag} is 'notransp'
##
## the result @code{A' * X}, if @var{flag} is 'transp'
##
## A typical case is @var{A} defined by @code{B ^ m}, 
## in which the result @code{A * X} can be computed without
## even forming explicitely @code{B ^ m} by 
##
## @code{Y = X;, for i = 1:m, Y = B * Y; endfor}
##
## The parameters @var{P1}, @var{P2}, ... are arguments of 
## @code{AFUN (FLAG, X, P1, P2, ...)}.
##
## The default value for @var{t} is 2. The algorithm requires
## matrix-matrix products with sizes @var{N}x@var{N} and
## @var{N}x@var{t}.
##
## The default initial matrix @var{X0} has the first column
## @code{ones(N,1) / N}
## and, if @var{t}>1, the remaining columns with random elements 
## @code{-1 / N}, @code{1 / N}.
##
## On output, @var{C} is the desired estimate, @var{V} and @var{W}
## vectors such that @code{W = A * V}, with
## @code{norm (W, 1)} = @code{C * norm (V, 1)}.
## @var{IT} contains in @code{IT(1)} the number of iterations
## (the maximum number is hardcoded to 5) and in  @code{IT(2)}
## the total number of products @code{A * X} or
## @code{A' * X} performed by the algorithm.
##
## Reference: @nospell{L.F. Shampine},
## @cite{A block algorithm for matrix 1-norm estimation, with and
## application to 1-norm pseudospectra}, SIAM J. Matrix Anal. Appl., 
## pp. 1185--1201, Vol 21, No. 4, 2000.
##
## @seealso{normest}
## @end deftypefn

function [est, v, w, k] = normest1(A, t, X, varargin)

  if ((nargin <= 1) || isempty (t))
    t = 2;
  end

  if isnumeric (A)
    ## A is a matrix
    Aisnum = true;
    n = length (A);
    if ((n <= 4) || (t == n))
    ## compute directly
      [est, idx] = max (sum (abs (A)));
      v = zeros(n, 1);
      v(idx) = 1;
      w = A(:, idx);
      k = [0, 1];
      return
    else
      realm = isreal (A);
    endif
  else
    Aisnum = false;
    n = A ('dim', [], varargin{:});
    realm = A ('real', [], varargin{:});
    Afun = @(x) A ('notransp', x, varargin{:});
    A1fun = @(x) A ('transp', x, varargin{:});
  endif

  if ((nargin <= 2) || isempty (X))
    ## save the seed, see normest
    ##    rs = rand ("state");
    ##    rand ("state", t);
    X = [ones(n, 1), sign(2 * rand(n, t - 1) - 1)]; 
    i = 2;
    imax = min (t, 2 ^ (n - 1));
    ##  There are at most 2^(n-1) unparallel columns, see later.
    while (i <= imax)
      if (any (abs (X(:,i)' * X(:,1:i-1)) == n))
        ## column i is parallel to a colum 1:i-1. Change it.
        X(:,i) = sign (2 * rand (n, 1) - 1);
      else
        i++;
      endif
    endwhile
    X = X / n;
  endif

  itmax = 5;
  ind_hist = zeros(n,1);
  estold = 0;
  ind = zeros (n, 1);
  S = zeros (n, t);
  k = [0; 0];
  conv = false;
  while ((~conv) && (k(1) < itmax))
    k(1)++;
    if (Aisnum)
      Y = A * X;
    else
      Y = Afun (X);
    endif
    k(2)++;
    [est, j] = max (sum (abs (Y)));
    if ((est > estold) || (k(1) == 2))
      ind_best = ind(j);    
      w = Y(:, j); ## there is an error in Algorithm 2.4
    endif
    if ((est <= estold) && (k(1) >= 2)) ## (1) of Algorithm 2.4
      est = estold;
      break; ## while
    endif
    estold = est;
    Sold = S;
    S = sign (Y);
    S(S==0) = 1;
    possiblebreak = false;
    if (realm)
      ## test parallel (only real case)
      if (all (any (abs (Sold' * S) == n))) ## (2) of Algorithm 2.4
        ## all columns of S parallel to a column of Sold, exit
        possiblebreak = true;
        conv = true;
      else
        if (t > 1)       
          ## at least two columns of S are not parallel
          i = 1;
          imax = min (t, 2 ^ (n - 1));
          while (i <= imax)
            ## The maximum number of parallel columns of length n with entries
            ## {-1,1} is 2^(n-1). Therefore, if the number of columns of S is
            ## greater than 2^(n-1), for sure some of them are parallel to some 
            ## columns of Sold. Don't even try to change them (i <= 2^(n-1)).
            ## Now, check if S(:,i) is parallel to any previous column of S
	    p = (any (abs (S(:,i)' * S(:,1:i-1)) == n));
            if (p || (any (abs (S(:,i)' * Sold) == n)))
              ## i-th column of S parallel to a previous or to a column of Sold:
              ## change it. 
	      S(:,i) = sign (2*rand (n, 1)-1);
	    else
	      i++;
	    endif
          endwhile
        endif
      endif
    endif
    if (~possiblebreak)
      if (Aisnum)
        Z = A' * S;
      else
        Z = A1fun (S); ## (3) of Algorithm 2.4
      endif
      k(2)++;
      h = max (abs (Z), [], 2);
      ind = (1:n)';
      if ((k(1) >= 2) && (max (h) == h(ind_best))) ## (4) of Algorithm 2.4
%        break; ## while
      endif
      [h, ind] = sort (h, "descend");
      if (t > 1)
        if (all (ind_hist(ind(1:t)))) ## (5) of Algorithm 2.4
          break ## while
        endif
        ind = ind(~ind_hist(ind));
        ## length(ind) could be less than t, especially if t is not << n.
        ## This is not considered in point (5) of Algorithm 2.4. 
        tmax = min (length (ind), t);
        ind = ind(1:tmax);
      else
        tmax = 1;
      endif
      X = zeros (n,tmax);
      for i = 1:tmax
        X(ind(i),i) = 1;
      endfor
      ind_hist(ind(1:tmax)) = 1;
    endif
  endwhile
  v = zeros (n, 1);    
  v(ind_best) = 1;
## restore state of random number generator
##  rand ("state", rs); 
  endfunction

%!function z = afun_A (flag, x, A, n)
%! switch flag
%! case {"dim"}
%!   z = n;
%! case {"real"}
%!   z = isreal (A);
%! case {"transp"}
%!   z = A' * x;
%! case {"notransp"}
%!   z = A * x;
%! endswitch
%!endfunction
%!function z = afun_A_P (flag, x, A, m)
%! switch flag
%! case {"dim"}
%!   z = length (A);
%! case {"real"}
%!   z = isreal (A);
%! case {"transp"}
%!   z = x; for i = 1:m, z = A' * z;, endfor
%! case {"notransp"}
%!   z = x; for i = 1:m, z = A * z;, endfor
%! endswitch
%!endfunction
%!test
%! ## simple test
%! A = reshape ((1:16)-8, 4, 4);
%! assert (normest1 (A), norm (A, 1), eps)
%!test
%! ## test t=1
%! A = rand (4); % for positive matrices always work
%! assert (normest1(A, 1), norm(A,1), 2 * eps)
%!test
%! ## test t=3
%! A = [-0.21825   0.16598   0.19388   0.75297;...
%!      -1.47732   0.78443  -1.04254   0.42240;...
%!       1.39857  -0.34046   2.28617   0.68089;...
%!       0.31205   1.50529  -0.75804  -1.22476];
%! X = [1,1,-1;1,1,1;1,1,-1;1,-1,-1]/3;
%! assert (normest1 (A, 3, X), norm (A, 1), 2 * eps)
%!test
%! ## test afun
%! A = rand (10);
%! n = length (A);
%! afun = @(flag, x) afun_A (flag, x, A, n);
%! assert (normest1 (afun), norm (A, 1), 2 * eps)
%!test
%! ## test afun with parameters
%! A = rand (10);
%! assert (normest1 (@afun_A_P, [], [], A, 3), norm (A ^ 3, 1), 2 * eps)
%!test
%! ## test output
%! A = reshape(1:16,4,4); 
%! [c, v, w, it] = normest1 (A);
%! assert(norm (w, 1), c * norm(v, 1), eps);
%! ## test output
%!test
%! A = rand (100);
%! A(A <= 1/3) = -1; 
%! A(A > 1/3 & A <= 2/3) = 0;
%! A(A > 2/3) = 1;
%! [c, v, w, it] = normest1 (A, 10);
%! assert(w, A * v, eps);
