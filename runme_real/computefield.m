function [STF,TVF] = computefield(u,Omega,sigma,rho)

%% PARAMETERS
%hsize        = 15;
hsize        = 80;
%scale_voting = [11:10:31];
scale_voting = [5:5:15];

% INLINE FUNCTIONS
blur      = @(f,sigma) imfilter(f,fspecial('gaussian',hsize,sigma),'replicate','same','conv');
nabla     = @(f)       grad_ij(f);
tensorize = @(u)       cat(4, cat(3,u(:,:,1).^2, u(:,:,1).*u(:,:,2)), cat(3,u(:,:,2).*u(:,:,1), u(:,:,2).^2));
T         = @(f,sigma) tensorize( nabla(blur(f,sigma)) );

normalize  = @(u) u./repmat(sqrt(sum(u.^2,3)), [1 1 2]);

delta      = @(S) (S(:,:,1,1)-S(:,:,2,2)).^2 + 4*S(:,:,1,2).^2;
eigenval1  = @(S) (S(:,:,1,1)+S(:,:,2,2)+sqrt(delta(S)))/2;
eigenval2  = @(S) (S(:,:,1,1)+S(:,:,2,2)-sqrt(delta(S)))/2;
eigenvec1  = @(S) normalize( cat(3, 2*S(:,:,1,2), S(:,:,2,2)-S(:,:,1,1)+sqrt(delta(S)) ) );
eigenvec2  = @(S) normalize( cat(3, -(S(:,:,2,2)-S(:,:,1,1)+sqrt(delta(S))), 2*S(:,:,1,2) ) );

recompose  = @(e1,e2,lambda1,lambda2) repmat(lambda1,[1 1 2 2]).*tensorize(e1) + repmat(lambda2,[1 1 2 2]).*tensorize(e2);

%% Compute Structure Tensor
ST = zeros(size(u,1),size(u,2),2,2);
TT = zeros(size(u,1),size(u,2),2,2);

for c = 1:size(u,3)
    
    Tloc           = T(u(:,:,c),sigma);
    STloc(:,:,1,1) = blur(Tloc(:,:,1),rho);
    STloc(:,:,1,2) = blur(Tloc(:,:,2),rho);
    STloc(:,:,2,1) = blur(Tloc(:,:,3),rho);
    STloc(:,:,2,2) = blur(Tloc(:,:,4),rho);
    
    
    TT = TT + STloc;
    ST = ST + STloc;
    
end

STF.e1 = eigenvec1(ST);
STF.e2 = eigenvec2(ST);
STF.l1 = eigenval1(ST);
STF.l2 = eigenval2(ST);

% get orientation of directions in (x,y)
ori_stf = -pi/2 + atan2(-STF.e1(:,:,1),STF.e1(:,:,2));
STF.e1  = cat(3,cos(ori_stf),sin(ori_stf));
STF.e2  = cat(3,-sin(ori_stf),cos(ori_stf));        

STF.E   = sqrt(STF.l1 + STF.l2);
STF.A   = (STF.l1 - STF.l2) ./ (STF.l1 + STF.l2);
STF.sal = (STF.l1 - STF.l2);

%% Compute Tensor Voting
TV = zeros(size(u,1),size(u,2),2,2);

for kk=1:size(u,3)
    
    % get saliency and orientation
    [sal_loc,ori_loc] = encode(blur(u(:,:,kk),sigma));
    %[sal_loc,ori_loc] = encode(u(:,:,kk));
    
    % zeroing the shadow edges information
    sal_loc = sal_loc.*Omega;
    ori_loc = ori_loc.*Omega + 2*pi.*(rand(size(ori_loc))-0.5).*(~Omega);  
    
    for ss=1:numel(scale_voting)
        
        [saliency,ballness,orientation] = vote(sal_loc,ori_loc,scale_voting(ss));
 
        l1 = saliency + ballness;
        l2 = ballness;
        % now orientation is in (x,y)
        e1 = cat(3,cos(orientation),sin(orientation));
        e2 = cat(3,-sin(orientation),cos(orientation));        

        TV = TV + recompose(e1,e2,l1,l2);
        
    end
    
end

TVF.e1 = eigenvec1(TV);
TVF.e2 = eigenvec2(TV);
TVF.l1 = eigenval1(TV);
TVF.l2 = eigenval2(TV);

% get orientation of directions in (x,y) coorinates
ori_tvf = -pi/2 + atan2(-TVF.e1(:,:,1),TVF.e1(:,:,2));
TVF.e1  = cat(3,cos(ori_tvf),sin(ori_tvf));
TVF.e2  = cat(3,-sin(ori_tvf),cos(ori_tvf));        

TVF.E   = sqrt(l1+l2);
TVF.A   = (l1-l2)./(l1+l2);
TVF.sal = (l1-l2);

%% RETURN THE GRADIENT
% TVF.z         = TVF.e2 is the gradient
% TVF.z_perp    = TVF.e1 is the main direction,
% STF.zold      = STF.e2 is the gradient
% STF.zold_perp = STF.e1 is the main direction,

% ATTENTION 
% TVF.e1(:,:,2) is the y-th component orientated along y
% TVF.e1(:,:,1) is the x-th component orientated along x
% STF.e1(:,:,2) is the y-th component orientated along y
% STF.e1(:,:,1) is the x-th component orientated along x
% atan2d(STV.e1e(:,:,2),STV.e1e(:,:,1)) returns the angle alpha in (x,y) coordinates
TVF.theta_z      = atan2d(TVF.e2(:,:,2),TVF.e2(:,:,1));
TVF.theta_z_perp = atan2d(TVF.e1(:,:,2),TVF.e1(:,:,1));
STF.theta_z      = atan2d(STF.e2(:,:,2),STF.e2(:,:,1));
STF.theta_z_perp = atan2d(STF.e1(:,:,2),STF.e1(:,:,1));

% project the angle in (i,j)
TVF.z       = cat(3,-sin(deg2rad(TVF.theta_z)),     cos(deg2rad(TVF.theta_z)));
TVF.z_perp  = cat(3,-sin(deg2rad(TVF.theta_z_perp)),cos(deg2rad(TVF.theta_z_perp)));
STF.z       = cat(3,-sin(deg2rad(STF.theta_z)),     cos(deg2rad(STF.theta_z)));
STF.z_perp  = cat(3,-sin(deg2rad(STF.theta_z_perp)),cos(deg2rad(STF.theta_z_perp)));

end

function DU = grad_xy(f)
% gradient in (x,y) coordinates

px = zeros(size(f));
py = zeros(size(f));
px(:,1:end-1) = f(:,2:end)-f(:,1:end-1);
py(2:end,:) = f(1:end-1,:)-f(2:end,:);

% in (x,y) coordinates
DU(:,:,1) = px;
DU(:,:,2) = py;


end

function DU = grad_ij(f)
% gradient in (i,j) coordinates

pi = zeros(size(f));
pj = zeros(size(f));
pi(1:end-1,:) = f(2:end,:)-f(1:end-1,:);
pj(:,1:end-1) = f(:,2:end)-f(:,1:end-1);

% in (i,j) coordinates
DU(:,:,1) = pi;
DU(:,:,2) = pj;


end