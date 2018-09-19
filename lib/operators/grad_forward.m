function [DU,D1,D2] = grad_forward(u)

% FORWARD DIFFERENCE
M = size(u,1);
N = size(u,2);

D1 = spdiags([-ones(M,1) ones(M,1)],[0 1],M,M);
D2 = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N);

% Boundary conditions
D1(M,:) = 0;
D2(N,:) = 0;


D1 = kron(speye(N),D1);
D2 = kron(D2,speye(M));
	
DU(:,:,1) = reshape(D1*u(:),M,N);
DU(:,:,2) = reshape(D2*u(:),M,N);

return