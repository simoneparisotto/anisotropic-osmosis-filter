function [DU,D1,D2] = grad_centered(u)

% Compute gradient as in
% - J. Weickert and H. Scharr. 
% "A Scheme for Coherence-Enhancing Diffusion Filtering with Optimized  
% Rotation  Invariance."
% http://dx.doi.org/10.1006/jvci.2001.049
% - P. Getreuer, 
% "Roussos-Maragos Tensor-Driven Diffusion for Image Interpolation"
% http://www.ipol.im/pub/art/2011/g_rmdi/article.pdf

M = size(u,1);
N = size(u,2);

%% D1 and D2 along dimensions
D1 = spdiags([-ones(M,1), ones(M,1)],[-1 1],M,M);
D2 = spdiags([-ones(N,1), ones(N,1)],[-1 1],N,N);

%% SUM MATRIX
% along i (u_{i+1,j} + u_{i-1,j}) 
m1x  = spdiags(ones(M,2),[-1,1],M,M);
% along j (u_{i,j+1} + u_{i,j-1})
m1y  = spdiags(ones(N,2),[-1,1],N,N);

%% Boundary conditions
D1([1 M],:) = 0;
D2([1 N],:) = 0;
%D1([1 2 3 M-2 M-1 M],:) = 0;
%D2([1 2 3 N-2 N-1 N],:) = 0;

m1x([1 M],:) = 0;
m1y([1 N],:) = 0;
%m1x([1 2 3 M-2 M-1 M],:) = 0;
%m1y([1 2 3 N-2 N-1 N],:) = 0;

%% build 2D operator
D1   = kron(speye(N),D1);
D2   = kron(D2,speye(M));
M1x  = kron(speye(N),m1x);
M1y  = kron(m1y,speye(M));

% D1 and D2
% D1 = (10*D1 + 3*M1y*D1)/32;
% D2 = (10*D2 + 3*M1x*D2)/32;
D1 = (10*D1 + 3*M1y*D1)/32;
D2 = (10*D2 + 3*M1x*D2)/32;

%% COMPUTE EXPLICIT DERIVATIVE (if necessary)
DU(:,:,1) = reshape(D1*u(:),M,N);
DU(:,:,2) = reshape(D2*u(:),M,N);

return