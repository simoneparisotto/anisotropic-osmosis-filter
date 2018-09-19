function extr_eigs = gersh (A, dense)
% function extr_eigs = gersh (A)
% function extr_eigs = gersh (A, 'sparse')
% function extr_eigs = gersh (A, 'full')
%
% Description: Compute extremal real and imaginary parts of eigenvalues.
%
% Compute the extremal largest real part (extr_eigs.LR), smallest real
% part (extr_eigs.SR), largest imaginary part (extr_eigs.LI), and
% smallest imaginary part (extr_eigs.SI) of the spectrum of A, by
% Gershgorin's disks. A is split into its Hermitian A_H and
% skew-Hermitian A_sH parts, before computation.
% Then, \sigma(A) \subseteq W(A) = W (A_H+A_sH) \subseteq W(A_H) +
% W(A_sH) = conv(\sigma(A_H)) + conv(\sigma(A_sH)) with
%
% conv(\sigma(A_H)) \subseteq [extr_eigs.SR,extr_eigs.LR]
% conv(\sigma(A_sH)) \subseteq i*[extr_eigs.SI,extr_eigs.LI]
%
% This file is part of the expleja project.
% Authors: Marco Caliari and Peter Kandolf
% Date: September 23, 2016
%
  if (((nargin == 1) && issparse (A)) || ...
      ((nargin == 2) && strcmp (dense, 'sparse')))
    center = 2 * diag (A);
  %% Hermitian part
    centers = real (center);
    radius = sum (abs (A + A'), 1)' - abs (centers);
    extr_eigs.SR = full (min (centers - radius)) / 2;
    if (norm (centers, inf) == 0)
      extr_eigs.LR = -extr_eigs.SR;
    else
      extr_eigs.LR = full (max (centers + radius)) / 2;
    end
    %% skew-Hermitian part
    centers = imag (center);
    radius = sum (abs (A - A'), 1)' - abs (centers);
    extr_eigs.SI = full (min (centers - radius)) / 2;
    if (norm (centers, inf) == 0)
      extr_eigs.LI = -extr_eigs.SI;
    else
      extr_eigs.LI = full (max (centers + radius)) / 2;
    end
  else
    radius = zeros (size (A, 2), 2);
    for i = 1:size (A, 2)
      radius(i, 1) = sum (abs (A(i,:)' + A(:,i)));
      radius(i, 2) = sum (abs (A(i,:)' - A(:,i)));
    end
    center = 2 * diag (A);
    centers = real (center);
    radius(:,1) = radius(:,1) - abs (centers);
    extr_eigs.SR = min (centers - radius(:,1)) / 2;
    if (norm (centers, inf) == 0)
      extr_eigs.LR = -extr_eigs.SR;
    else
      extr_eigs.LR = max (centers + radius(:,1)) / 2;
    end
    centers = imag (center);
    radius(:,2) = radius(:,2) - abs (centers);
    extr_eigs.SI = min (centers - radius(:,2)) / 2;
    if (norm (centers, inf) == 0)
      extr_eigs.LI = -extr_eigs.SI;
    else
      extr_eigs.LI = max (centers + radius(:,2)) / 2;
    end
  end
end
%!test
%! A = toeplitz([-2,1,0,0])+toeplitz([0,1,0,0],[0,-1,0,0]);
%! extr_eigs = gersh (A);
%! ref.SR = -4;
%! ref.LR = 0;
%! ref.LI = 2;
%! ref.SI = -2;
%! assert(extr_eigs, ref)
%!test
%! A = toeplitz([-2,1,0,0])+toeplitz([0,1,0,0],[0,-1,0,0]);
%! extr_eigs = gersh (A, 'sparse');
%! ref.SR = -4;
%! ref.LR = 0;
%! ref.LI = 2;
%! ref.SI = -2;
%! assert(extr_eigs, ref)
%!test
%! A = sparse (toeplitz([-2,1,0,0])+toeplitz([0,1,0,0],[0,-1,0,0]));
%! extr_eigs = gersh (A);
%! ref.SR = -4;
%! ref.LR = 0;
%! ref.LI = 2;
%! ref.SI = -2;
%! assert(extr_eigs, ref)
%!test
%! A = sparse (toeplitz([-2,1,0,0])+toeplitz([0,1,0,0],[0,-1,0,0]));
%! extr_eigs = gersh (A, 'full');
%! ref.SR = -4;
%! ref.LR = 0;
%! ref.LI = 2;
%! ref.SI = -2;
%! assert(extr_eigs, ref)
%!test
%! A = toeplitz([-2,1,0,0])+toeplitz([0,1,0,0],[0,-1,0,0]);
%! A = 1i * A;
%! extr_eigs = gersh (A);
%! ref.SR = -2;
%! ref.LR = 2;
%! ref.LI = 0;
%! ref.SI = -4;
%! assert(extr_eigs, ref)
%!test
%! A = toeplitz([-2,1,0,0])+toeplitz([0,1,0,0],[0,-1,0,0]);
%! extr_eigs = gersh (sparse (A));
%! ref.SR = -4;
%! ref.LR = 0;
%! ref.LI = 2;
%! ref.SI = -2;
%! assert(extr_eigs, ref)
%!test % real output
%! A = 2*sprand(100,100,0.1)-1+1i*(2*sprand(100,100,0.1));
%! extr_eigs = gersh (A);
%! assert(isreal(extr_eigs.LR) && isreal(extr_eigs.SR) && ...
%! isreal(extr_eigs.SI) && isreal(extr_eigs.LI),true)
%!test % opposite LI and SI for real input
%! A = 2*sprand(100,100,0.1)-1;
%! extr_eigs = gersh (A);
%! assert(extr_eigs.LI+extr_eigs.SI,0)
%!test % opposite LR and SR for pure imaginary input
%! A = 1i*(2*sprand(100,100,0.1)-1);
%! extr_eigs = gersh (A);
%! assert(extr_eigs.LR+extr_eigs.SR,0)
%!test % LI == SI == 0 for hermitian input
%! A = 2*sprand(100,100,0.2)-1+1i*(2*sprand(100,100,0.2)-1);
%! A = A*A';
%! extr_eigs = gersh (A);
%! assert((extr_eigs.LI == 0) && (extr_eigs.SI == 0),true)
