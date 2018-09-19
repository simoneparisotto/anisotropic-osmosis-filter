%function [P3, P3X_old, P3Y_old, P3_old] = directional_osmosis_discretization(umat,umask,WW)
function P3 = coherence_osmosis_discretization(umat,umask,WW)
% Input:
%   umat  = shadowed image
%   umask = shadow boundary indicator function
%   WW    = anisotropic diffusion tensor m x n x 4

[mx, my, c] = size(umat);

x  = linspace(1,mx,mx)';
y  = linspace(1,my,my)';
hx = (max(x)-min(x))/(mx-1);
hy = (max(y)-min(y))/(my-1);

%% DIRECTIONAL WEIGHTS
%   [1 0; 0 1] outside the shadow border
%   b1 v \otimes v + b2 vpepr \otimes vperp on the shadow border
W{1} = WW(:,:,1);
W{2} = WW(:,:,2);
W{3} = WW(:,:,3);
W{4} = WW(:,:,4);

%% FIRST DERIVATIVES
[~,D1,D2]   = grad_centered(ones(mx,my));

%% GENERALIZED OSMOSIS FILTER MASK
mask = double(umask(:)==1);   % SHADOW REMOVAL

%% ANISOTROPIC DIFFUSION
% diffusion = div(W \grad u)
W1 = spdiags(W{1}(:),0,size(D1,1),size(D1,2))*D1;
W2 = spdiags(W{2}(:),0,size(D2,1),size(D2,2))*D2;
W3 = spdiags(W{3}(:),0,size(D1,1),size(D1,2))*D1;
W4 = spdiags(W{4}(:),0,size(D2,1),size(D2,2))*D2;

diffusion = -(D1.'*(W1+W2) + D2.'*(W3+W4));

%% DRIFT
% drift = \grad(log(f)) for umask = 1
% drift = 0             for umask = 0
%d1 = mask.*D1*log(umat(:));
%d2 = mask.*D2*log(umat(:));
d1 = mask.*(D1*umat(:))./umat(:);
d2 = mask.*(D2*umat(:))./umat(:);


% compute WD = [WD1,WD2]
WD11 = reshape( W{1}(:).*d1 ,size(umat) );
WD12 = reshape( W{2}(:).*d2 ,size(umat) );
WD21 = reshape( W{3}(:).*d1 ,size(umat) );
WD22 = reshape( W{4}(:).*d2 ,size(umat) );
Wd1 = WD11 + WD12;
Wd2 = WD21 + WD22;
WD1 = spdiags(Wd1(:),0,size(D1,1),size(D1,2));
WD2 = spdiags(Wd2(:),0,size(D2,1),size(D2,2));

% drift = div(Wd u)
drift = -(D1.'*WD1 + D2.'*WD2);

%% ARTIFICIAL DIFFUSION
% based on Quarteroni, Numerical Models for Differential Problems, p.162
% diffusion_streamline = -Q*h*div( (B\cdot\grad u) B );
% our model reads in space (Th. 4.5 part 2):
% div(W(\grad u - du)) = W\cdot D^2u + (div(W)-Wd)\cdot\grad u - div(Wd)u
% so B = [B1,B2] = [div(W)-Wd]
flag_artificial = 0;
if flag_artificial
    % pixel size
    h = 1;
    
    % div(W), W is m x n x 4
    divW1 = reshape( -(D1.'*W{1}(:)+D2.'*W{2}(:)) ,size(umat));
    divW2 = reshape( -(D1.'*W{3}(:)+D2.'*W{4}(:)) ,size(umat));
    
    % B = div(W)-Wd
    B1 = divW1-Wd1;
    B2 = divW2-Wd2;
    
    % coefficients for (B\cdot\grad u) B
    B11 = B1.*B1;
    B12 = B1.*B2;
    B21 = B2.*B1;
    B22 = B2.*B2;
    BB11 = spdiags(B11(:),0,size(D1,1),size(D1,2));
    BB12 = spdiags(B12(:),0,size(D1,1),size(D1,2));
    BB21 = spdiags(B21(:),0,size(D2,1),size(D2,2));
    BB22 = spdiags(B22(:),0,size(D2,1),size(D2,2));
    
    % A = (B\cdot\grad u) B
    A1 = BB11*D1 ;%+ BB12*D2;
    %A2 = BB21*D1 +
    A2 = BB22*D2;
    
    % div( (B\cdot\grad u) B )
    divB = -(D1.'*A1 + D2.'*A2);
    
    % Q = epsilon(x,y) (for now).. it should be
    % Q = 1./abs(B1+B2);
    AA = 10000*sqrt(sum(cat(3,reshape(A1*umat(:),size(umat)),reshape(A1*umat(:),size(umat))).^2,3));
    %AA(AA<0)=0;
    %AA(isinf(AA) | isnan(AA)) = 0;
    % but I get +\infty
    Q = AA.*ones(size(umat)); % for =0 NO ARTIFICIAL DIFFUSION
    Q = spdiags(Q(:),0,size(D1,1),size(D1,2));
    divB(divB<0)=0;
    diffusion_art = Q*h*divB;
end

% %% CLASSIC
% outshadow_1_pre    = padarray(umask>0,[1,0,0],'replicate','pre');
% outshadow_2_pre    = padarray(umask>0,[0,1,0],'replicate','pre');
% 
% [~,D1classic,D2classic] = grad_forward(ones(mx,my));
% 
% % average upper (u_{i+1,j} + u_{ij})/2  and lower (u_{ij}+u_{i-1,j})/2
% m1xup  = spdiags(ones(mx,2)/2,[0,1],mx,mx);
% m1xlow = spdiags(ones(mx,2)/2,[-1,0],mx,mx);
% % average upper (u_{i,j+1} + u_{ij})/2  and lower (u_{ij}+u_{i,j-1})/2
% m1yup  = spdiags(ones(my,2)/2,[0,1],my,my);
% m1ylow = spdiags(ones(my,2)/2,[-1,0],my,my);
% 
% M1xup  = kron(speye(my),m1xup);
% M1xlow = kron(speye(my),m1xlow);
% M1yup  = kron(m1yup,speye(mx));
% M1ylow = kron(m1ylow,speye(mx));
% 
% % LAPLACIAN
% Dxx_old = -D1classic.'*D1classic;
% Dyy_old = -D2classic.'*D2classic;
% 
% % STANDARD DRIFT VECTOR FIELD d
% d1ij = zeros(mx+1,my,c);
% d2ij = zeros(mx,my+1,c);
% 
% d1ij(2:mx,:,:) = diff(umat,1,1)./(umat(2:mx,:,:)+umat(1:mx-1,:,:))*2/hx;
% d2ij(:,2:my,:) = diff(umat,1,2)./(umat(:,2:my,:)+umat(:,1:my-1,:))*2/hy;
% 
% % APPLY SHADOW MASK: DRIFT VECTOR = 0 across shadow boundary i.e. umask=0;
% d1ij(~outshadow_1_pre) = 0;
% d2ij(~outshadow_2_pre) = 0;

for k = (1:c)
    
    % GENERALIZED OSMOSIS FILTER
    P3{k} = diffusion - drift;
    if flag_artificial
        P3{k} = P3{k} + diffusion_art;
    end
    
%     %  STANDARD OSMOSIS FILTER
%     P3X_old{k} = Dxx_old - 1/hx*( ...
%         spdiags(reshape(d1ij(2:mx+1,:,k),mx*my,1),0,mx*my,mx*my)*M1xup - ...
%         spdiags(reshape(d1ij(1:mx,:,k),mx*my,1),0,mx*my,mx*my)*M1xlow);
%     
%     P3Y_old{k} = Dyy_old - 1/hy*( ...
%         spdiags(reshape(d2ij(:,2:my+1,k),mx*my,1),0,mx*my,mx*my)*M1yup - ...
%         spdiags(reshape(d2ij(:,1:my,k),mx*my,1),0,mx*my,mx*my)*M1ylow);
%     
%     
%     P3_old{k}  =  P3X_old{k} + P3Y_old{k};
   
    
end


end