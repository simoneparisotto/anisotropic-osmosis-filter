% Demo examples for expleja
% See %!demo section at the end of expleja.m
%
% Authors: Marco Caliari, Peter Kandolf
% Date: September 29, 2016
fprintf('1) Normal functionality of EXPLEJA:\n');
n = 10;
A = -gallery ('poisson', n);
v = linspace (-1, 1, n ^ 2)';
t = 1;
y = expleja (t, A ,v);
fprintf('\tRelative differences between EXPM and EXPLEJA.\n')
fprintf('\tShould be of order %9.2e and has error %9.2e\n', eps / 2, ...
        norm (y - expm( t * A ) * v, 1) / norm (y, 1))
fprintf('2) With multiple vectors v\n');
n = 10;
A = -gallery ('poisson', n);
v = linspace (-1, 1, n ^ 2)';
v1 = exp (-10 * v .^ 2);
t = 1;
y = expleja (t, A, [v, v1]);
fprintf('\tRelative differences between EXPM and EXPLEJA.\n')
fprintf('\tShould be of order %9.2e and has error %9.2e\n', eps / 2, ...
        norm (y - expm (t * A) * [v, v1], 1) / norm (y, 1))
fprintf('3) Precompute the interpolation parameters\n');
n = 10;
A = -gallery ('poisson', n);
v = linspace (-1, 1, n ^ 2)';
t = 1;
y1 = v;
mult = 10;
tic
for i = 1:mult
  y1 = expleja (t, A, y1);
end
time0 = toc;
% Standard configuration for expleja and param
extreigs = gersh (A);     % Gerschgorin estimates
tol = [0, 2 ^ -53, 1, 1]; % tolerance
p = 5;                    % maximal order of hump reduction
max_points = 100;         % maximal number of interpolation points
shift = 1;                % allow a matrix shift
tic
[param, ~] = select_interp_para (t, A, extreigs, tol, max_points, p, shift);
time1 = toc;
y2=v;
tic
for i = 1:mult
  y2 = expleja (t, A, y2, tol, p, param);
end
time2 = toc;
fprintf('\tTiming comparison in seconds with constant stepsize.\n')
fprintf('\tTime for %d steps w/o precomputing %9.2e\n', mult, time0)
fprintf('\tTime for %d steps with precomputing %9.2e\n', mult, time2)
fprintf('\tTime for the precomputing stage %9.2e\n', time1)
fprintf('\tTotal speed-up %9.2e\n', time0/(time1+time2))
fprintf('4) Simulate integrator with variable time step size\n');
A = -10 * gallery ('triw', 20, 4);
v = ones (length (A), 1);
t = 1;
y1 = v;
T = [1,1.25,1.1,1.5,1,1,1.25,1.1,1.5,1,1,1.25,1.1,1.5,1,1,1.25,1.1,1.5,1];
tic
for t = T
  y1 = expleja (t, A, y1);
end
time1 = toc;
% Standard configuration for expleja and param
extreigs = gersh (A);     % Gerschgorin estimates
tol = [0, 2 ^ -53, 1, 1]; % tolerance
p = 5;                    % maximal order of hump reduction
max_points = 100;         % maximal number of interpolation points
shift = 1;                % allow a matrix shift
y2 = v;
tic
[param, ~] = select_interp_para(T(1), A, extreigs, tol, max_points, p,...
                               shift);
y2 = expleja (T(1), A, y2, tol, p, param);
shift = 0; % Do not reshift A in the param update
           % (ATTENTION keep p constant)
for t = T(2:end)
  [param, ~] = select_interp_para (t, A, extreigs, tol, max_points, p,...
                                   shift,param);
  y2 = expleja (t, A, y2, tol, p, param);
end
time2 = toc;
fprintf('\tTiming comparison in seconds with variable stepsize.\n')
fprintf('\tTime for multiple steps w/o precomputing %9.2e and error %9.2e\n',time1, norm(y1-expm(sum(T)*A)*v,1)/norm(y1,1))
fprintf('\tTime for multiple steps with precomputing %9.2e and error %9.2e\n',time2, norm(y2-expm(sum(T)*A)*v,1)/norm(y2,1))
fprintf('\tTotal speed-up %9.2e\n',time1/time2)
