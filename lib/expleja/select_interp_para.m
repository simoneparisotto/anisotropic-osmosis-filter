function [param, A] =...
        select_interp_para (H, A, extreigs, tol, m_max, p, shift, paramin)
% Function file: [param, A] =...
%                 select_interp_para (H, A, extreigs, tol, m_max, ...
%                 p, shift, paramin)
%
% Description: Compute parameters for matrix exponential interpolation.
%
% Computes the parameters for the interpolation, input is the time step H
%   the matrix A as well as the spectral estimate EXTREIGS.
%   TOL is used to select the correct parameters. M_MAX specifies the
%   maximal degree of interpolation (should be <=100) and P the maximal
%   power of A used to extimate the spectral radius. The boolean variable
%   SHIFT indicates if the input matrix A should be shifted or not. An
%   additional paramter PARAMIN can be given. If so the parameters are
%   updated according to H, TOL and P.
%
% The output structure param contains the number of substeps NSTEPS,
%   the interpolation interval GAMMA2, the scaled interpolation points XI,
%   the divided differences DD, the shift MU, the amount of matrix-vector
%   products required for the norm computation C, the suggested function to
%   perform the newton interpolation NEWT, the estimated degree
%   of interpolation M as well as the computed estimates for ||A^p||^(1/p)
%   as dest. The (shifted) matrix A is returned separately.
%
% This function can be called independently from EXPLEJA to precompute
%   PARAM. This might have computational advantages if A does not change
%   during successive computations.
% If an update of these parameters is required for different H parma can be
%   a input parameter.
%
% This file is part of the expleja project.
% Authors: Marco Caliari and Peter Kandolf
% Date: September 23, 2016
  h = abs (H);
  n = size (A, 2);
  %% Check if param needs to be updated or newly computed
  if nargin>=8 && ~isempty(paramin) && isstruct(paramin)
    param=paramin;  
    if (paramin.h==H) && sum(tol==paramin.tol)==4
      param.c=0;
      if shift
        A=A-param.mu*speye(n,n);
      end
      return
    end
  else
    paramin=[];
  end
  %% Shift for the matrix
  mu=(extreigs.SR+extreigs.LR)/2+1i*(extreigs.SI+extreigs.LI)/2;
  %% Select correct method depending on the spectral estimate
  %% (real vs. complex conjugate)
  isrealint = true;
  if abs (extreigs.LR - extreigs.SR)  >= abs (extreigs.LI - extreigs.SI)
    newt = @newton;
  else
    isrealint = false;
    newt = @newtons;
  end
  %% Rename A to save storage
  if shift
    A = A - mu * speye (n, n); % if A is diagonal constant, then normA == 0
  end
  %% Select tolerance and corresponding file with stored data
  sampletol = [2 ^ -10, 2 ^ -24, 2 ^ -53];
  tol_strings = {'half', 'single','double'};
  t = min (tol(tol(1:2) ~= 0));
  if isempty(t) || t < 2 ^ -53, t = 2 ^ -53; end
  %% Load file information
  index = find (sampletol <= t, 1, 'first');
  if isempty (index), index = length (sampletol); end
  if isrealint
    filename = sprintf('data_leja_%s.mat', tol_strings{index});
  else
    filename = sprintf('data_lejas_%s.mat', tol_strings{index});
    m_max = floor (m_max / 2);
  end
  file = load (filename);
  %% Compute estimate of \|A^i\|^(1/i) for a hump reduction, if beneficial
  c=0;
  if ~isempty(paramin) && paramin.tol(4)==tol(4) && ...
            length(paramin.dest)==p
        dest=paramin.dest;
        normA=dest(1);
        reduction = paramin.reduction;
  else
    if ~isempty(paramin) && paramin.tol(4)==tol(4)
        normA=paramin.dest(1);
        c=0;
    else
        [normA,c] = normAmp(1,tol(4));
    end
    dest = []; 
    reduction = false;
    if (p > 0 && normA ~= 0)
        dest = inf (p,1); dest(1) = normA; j = 2;
        if ~isempty(paramin) && paramin.tol(4)==tol(4) && ...
                length(paramin.dest)>=j
            dest(j)=paramin.dest(j);
        else
            [dest(j),k]=normAmp(j,tol(4)); c=c+k;
        end
        while (j + 1 <= p && dest(j - 1) / dest(j) > 1.01)
            reduction = true;
            j = j + 1; 
            if ~isempty(paramin) && paramin.tol(4)==tol(4) && ...
                 length(paramin.dest)>=j
                dest(j)=paramin.dest(j);
            else
                 [dest(j),k]=normAmp(j,tol(4)); c=c+k;
            end
         end
    end
  end
  theta = file.theta;
  %% Substep estimate for circles
  mm = min (m_max, length (theta));
  [nsteps, m] = min (ceil ((h * normA) ./ theta(1:mm)));
  %% make sure nsteps >=1 and m >= 2
  nsteps = max (nsteps, 1);
  k = m;
  %% Substep estimates for ellipses
  e = file.ell_eps;
  haxis = file.haxis;
  a = h * (extreigs.LR-extreigs.SR)/2;
  b = h * (extreigs.LI-extreigs.SI)/2;
  S = zeros (length (haxis), 3);
  if (reduction)
    for j=1:length(haxis)
      l=ceil(sqrt([(a+e)^2,(b+e)^2]*haxis{j}(:,1)));
      S(j,:)=[l,1,l*j];
    end
  else

    for j=1:length(haxis)
      [l,o]=min(ceil(sqrt([(a+e)^2,(b+e)^2]*(haxis{j}))));
      S(j,:)=[l,o,l*j];
    end
  end
  [~,j]=min(S(:,3));
  nsteps_ellipse=S(j,1);
  m_ellipse=j;
  if S(j,2)>1
    k_ellipse=S(j,2)+j-2;
  else
    tmp=1./sqrt(haxis{j}(1));
    k_ellipse=find(theta>=tmp(1),1,'first');
  end
  %% Select minimum of the two estimates
  if (nsteps * m > nsteps_ellipse * m_ellipse)
    nsteps = nsteps_ellipse;
    m = m_ellipse;
    k = k_ellipse;
  end
  if ~isrealint
    m = 2 * m;
  end
  %% Perform reduction if indicated
  if reduction
    gamma2 = h / nsteps * min(dest);
    k = find (theta >= gamma2, 1, 'first');
  end
  %Select correct parameters
  gamma2 = theta(k);
  dd = file.dd(:,k);
  dd = dd(1:m+1);
  xi = file.xi(1:m+1) * (gamma2 / 2);
  param = struct('nsteps', nsteps, 'gamma2', gamma2, 'xi', xi, 'dd', dd,...
               'mu', mu, 'c', c, 'newt', newt, 'm', m, 'dest', dest, ...
               'h', H, 'tol', tol,'reduction',reduction);
  %% nested function for handlying the norms
  function [c, mv] = normAmp (p, norm_type)
    switch norm_type
      case 1
        if  p == 1, c = norm (A, 1); mv = 0;
        else [c, ~, ~, tmv] = normest1(@(flag,x)normfun(flag,x,A,p),1);
          c = c ^ (1 / p); mv = tmv(2) * p;
        end
      case 2
        if (length(A) > 50)
           % MARCO why don't you do the same for the 1-norm? (that is,
           % compute explicitely norm (full(A)^m, 1) if length(A) < 50
          [c, mv] = normest2(@(flag,x)normfun(flag,x,A,p),[],8);
          c = c ^ (1 / p); mv = mv * p;
        else
          c = norm (full(A) ^ p, 2) ^ ( 1 / p); mv = 0;
        end
      case inf
        if p == 1, c = norm (A, inf); mv = 0;
        else
          [c,~,~,tmv] = normest1(@(flag,x)normfun(flag,x,A',p),1);
          c = c ^ (1 / p); mv = tmv(2) * p;
        end
    end
  end
end
%% subfunction
function x = normfun (flag, x , A ,m)
  switch flag
    case 'dim'
      x = size (A, 1);
    case 'real'
      x = isreal (A);
    case 'notransp'
      for i = 1:m, x = A * x; end
    case 'transp'
      for i = 1:m, x = A' * x; end
  end
end
