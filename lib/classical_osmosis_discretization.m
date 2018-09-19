function [P3X_old, P3Y_old, P3_old] = classical_osmosis_discretization(umat,umask)
% Input:
%   umat  = shadowed image
%   umask = shadow boundary indicator function

[mx, my, c] = size(umat);

x  = linspace(1,mx,mx)';
y  = linspace(1,my,my)';
hx = (max(x)-min(x))/(mx-1);
hy = (max(y)-min(y))/(my-1);

%% CLASSIC
outshadow_1_pre    = padarray(umask>0,[1,0,0],'replicate','pre');
outshadow_2_pre    = padarray(umask>0,[0,1,0],'replicate','pre');

[~,D1classic,D2classic] = grad_forward(ones(mx,my));

% average upper (u_{i+1,j} + u_{ij})/2  and lower (u_{ij}+u_{i-1,j})/2
m1xup  = spdiags(ones(mx,2)/2,[0,1],mx,mx);
m1xlow = spdiags(ones(mx,2)/2,[-1,0],mx,mx);
% average upper (u_{i,j+1} + u_{ij})/2  and lower (u_{ij}+u_{i,j-1})/2
m1yup  = spdiags(ones(my,2)/2,[0,1],my,my);
m1ylow = spdiags(ones(my,2)/2,[-1,0],my,my);

M1xup  = kron(speye(my),m1xup);
M1xlow = kron(speye(my),m1xlow);
M1yup  = kron(m1yup,speye(mx));
M1ylow = kron(m1ylow,speye(mx));

% LAPLACIAN
Dxx_old = -D1classic.'*D1classic;
Dyy_old = -D2classic.'*D2classic;

% STANDARD DRIFT VECTOR FIELD d
d1ij = zeros(mx+1,my,c); 
d2ij = zeros(mx,my+1,c);

d1ij(2:mx,:,:) = diff(umat,1,1)./(umat(2:mx,:,:)+umat(1:mx-1,:,:))*2/hx;
d2ij(:,2:my,:) = diff(umat,1,2)./(umat(:,2:my,:)+umat(:,1:my-1,:))*2/hy;

% APPLY SHADOW MASK: DRIFT VECTOR = 0 across shadow boundary i.e. umask=0;
d1ij(~outshadow_1_pre) = 0;
d2ij(~outshadow_2_pre) = 0;

for k = (1:c)
    
    %  STANDARD OSMOSIS FILTER                   
    P3X_old{k} = Dxx_old - 1/hx*( ...
        spdiags(reshape(d1ij(2:mx+1,:,k),mx*my,1),0,mx*my,mx*my)*M1xup - ...
        spdiags(reshape(d1ij(1:mx,:,k),mx*my,1),0,mx*my,mx*my)*M1xlow);
    
    P3Y_old{k} = Dyy_old - 1/hy*( ...
        spdiags(reshape(d2ij(:,2:my+1,k),mx*my,1),0,mx*my,mx*my)*M1yup - ...
        spdiags(reshape(d2ij(:,1:my,k),mx*my,1),0,mx*my,mx*my)*M1ylow);
   
    P3_old{k}  =  P3X_old{k} + P3Y_old{k};                
   
end


end