function d = exptaylordd (xi)
%
% Function file: d = exptaylordd(xi)
%
% Description: Compute divided differences for exp in a stable way.
%
% divided differences at xi points of exp
% it is equivalent (but much more accurate) to
%
% d = zeros(1,m);
% for i = 1:m
%   d(i) = exp (xi(i));
%   for j = 1:i - 1
%     d(i) = (d(i) - d(j)) / (xi(i) - xi(j));
%   end
% end
%
m = length (xi);
s = max(1,ceil(max(abs(xi))/0.7)); % McCurdy, Ng and Parlett 1984
F = toeplitz(1 ./ s .^ (0:m-1)./[1,1,cumprod(2:m-1)]);
for l = 1:16
% Matlab code
  for j = 1:m-1
    F(j, j) = xi(j) * F(j, j) / l / s;
    for i = j+1:m
      F(i, j) = (xi(i) * F(i, j) + F(i - 1,j)) / (l + i - j) / s;
      F(j, i) = F(j, i) + F(i, j);
    end
%    F(j, j + 1:m) = F(j, j + 1:m) + F(j + 1:m, j).';
  end
% end of Matlab code
% GNU Octave code
% Since GNU Octave has no jit, the following might be faster
%  A = spdiags([-1 ./ (l+1:l+m)' / s,ones(m,1)],[-1,0],m,m);
%  A(1,1) = 2;
%  A(2,1) = 2*A(2,1);
%  for j = 1:m-1
%    b = xi(j:m)(:) ./ (l:l+m-j)' / s .* F(j:m,j);
%    F(j:m,j) = A(1:m+1-j,1:m+1-j)\b;
%  end
%  F = F + tril (F).';
% end of GNU Octave code
end
F(1:m+1:m^2) = exp (xi / s);
F = triu(F);
d = F(1, :);
for i = 1:s-1
    d = d * F;
end
d=reshape(d,size(xi));
%!test
%! m = 3;
%! xi = 4*[-2,2,0];
%! d = zeros(1,m);
%! for i = 1:m
%!   d(i) = exp (xi(i));
%!   for j = 1:i - 1
%!     d(i) = (d(i) - d(j)) / (xi(i) - xi(j));
%!   end
%! end
%! d1 = exptaylordd (xi);
%! assert (max (abs (d - d1)) < 1e-14 * max (abs (d1)))
%!test
%! m = 3;
%! xi = 0*[-2,2,0];
%! d = exptaylordd (xi);
%! assert (max (abs (d - [1,1,1/2])) < 1e-14 * max (abs (d)))
