function [P3] = classical_osmosis_discretization_W(umat,umask,WW)
% Input:
%   umat  = shadowed image
%   umask = shadow boundary indicator function
%   WW : directional weights on the shadow boundaries

[mx, my, c] = size(umat);

x  = linspace(1,mx,mx)';
y  = linspace(1,my,my)';
hx = (max(x)-min(x))/(mx-1);
hy = (max(y)-min(y))/(my-1);

% outshadow_1_pre    = padarray(umask>0,[1,0,0],'replicate','pre');
% outshadow_2_pre    = padarray(umask>0,[0,1,0],'replicate','pre');
% outshadow_1_pre = outshadow_1_pre(1:end-1,:,:);
% outshadow_2_pre = outshadow_2_pre(:,1:end-1,:);

%% DIRECTIONAL WEIGHTS
%   [1 0; 0 1] outside the shadow border
%   b1 v \otimes v + b2 vpepr \otimes vperp on the shadow border
WLAP{1} = WW(:,:,1);
WLAP{2} = WW(:,:,2);
WLAP{3} = WW(:,:,3);
WLAP{4} = WW(:,:,4);

%% FIRST DERIVATIVES
[~,D1diff,D2diff]   = grad_forward(ones(mx,my));
[~,D1drift,D2drift] = grad_forward(ones(mx,my));

%% GENERALIZED OSMOSIS FILTER
mask = double(umask(:)==1);  % SHADOW REMOVAL

% diff = div(W grad u)
W12d = spdiags(WLAP{1}(:),0,size(D1diff,1),size(D1diff,2))*D1diff + spdiags(WLAP{2}(:),0,size(D2diff,1),size(D2diff,2))*D2diff;
W34d = spdiags(WLAP{3}(:),0,size(D1diff,1),size(D1diff,2))*D1diff + spdiags(WLAP{4}(:),0,size(D2diff,1),size(D2diff,2))*D2diff;

diffusion = - D1diff.'*W12d - D2diff.'*W34d;

% drift = div(W grad(log(f)))
d1 = (D1drift*umat(:))./umat(:);
d2 = (D2drift*umat(:))./umat(:);

dr1 = spdiags(mask.*(WLAP{1}(:).*d1 + WLAP{2}(:).*d2),0,size(D1diff,1),size(D1diff,2));
dr2 = spdiags(mask.*(WLAP{3}(:).*d1 + WLAP{4}(:).*d2),0,size(D2diff,1),size(D2diff,2));

%dr1 = spdiags(outshadow_1_pre(:).*WLAP{1}(:).*d1 + outshadow_1_pre(:).*WLAP{2}(:).*d2,0,size(D1diff,1),size(D1diff,2));
%dr2 = spdiags(outshadow_2_pre(:).*WLAP{3}(:).*d1 + outshadow_2_pre(:).*WLAP{4}(:).*d2,0,size(D2diff,1),size(D2diff,2));


drift = -D1drift.'*dr1 - D2drift.'*dr2;

% ARTIFICIAL DIFFUSION
% first
% b1_ART = cat(3,reshape(D1diff*WLAP{1}(:) + D2diff*WLAP{2}(:),size(umat)),...
%                reshape(D1diff*WLAP{3}(:) + D2diff*WLAP{4}(:),size(umat)));
%           
% b2_ART = cat(3,reshape(mask.*(WLAP{1}(:).*d1 + WLAP{2}(:).*d2),size(umat)),...
%                reshape(mask.*(WLAP{3}(:).*d1 + WLAP{4}(:).*d2),size(umat)));
%            
% b = b1_ART + b2_ART;
% Q = 1./sqrt(b(:,:,1).^2 + b(:,:,2).^2);

% second
% A11 = -D1diff*D1diff;
% A12 = -D2diff*D1diff;
% A21 = -D1diff*D2diff;
% A22 = -D2diff*D2diff;
% Q = 1./sqrt(sum(WW(:,:,1).^2+WW(:,:,2).^2+WW(:,:,3).^2+WW(:,:,4).^2,3);
% 
% A = A11+A12+A21+A22

% % third
% bG = spdiags(reshape(V(:,:,1),[],1),0,size(D1diff,1),size(D1diff,2))*D1diff + spdiags(reshape(V(:,:,2),[],1),0,size(D1diff,1),size(D1diff,2))*D2diff;
% V1bG = spdiags(reshape(V(:,:,1),[],1),0,size(D1diff,1),size(D1diff,2)).*bG;
% V2bG = spdiags(reshape(V(:,:,2),[],1),0,size(D1diff,1),size(D1diff,2)).*bG;
% 
% Q = 1./abs(sum(V,3));
% Q(isinf(Q)) = 0;
% Q = spdiags(reshape(Q,[],1),0,size(D1diff,1),size(D1diff,2));
% A = (-Q*(-D1diff*V1bG -D2diff*V2bG));

% % third
% WD = cat(3,WLAP{1}.*reshape(d1,size(umat)) + WLAP{2}.*reshape(d2,size(umat)),...
%            WLAP{3}.*reshape(d1,size(umat)) + WLAP{4}.*reshape(d2,size(umat)));
%        
% b1 = spdiags(reshape(WD(:,:,1),[],1),0,size(D1diff,1),size(D1diff,2))*D1diff;
% b2 = spdiags(reshape(WD(:,:,2),[],1),0,size(D1diff,1),size(D1diff,2))*D2diff;
% 
% divb = -(D1diff*b1 + D2diff*b2);
%        
% Q = 1./abs(WD(:,:,1)+WD(:,:,2));
% Q(Q>1) = 0;
% Q = spdiags(reshape(Q,[],1),0,size(D1diff,1),size(D1diff,2));
% A = -Q.*divb;

for k = (1:c)
    
     % GENERALIZED OSMOSIS FILTER
    P3{k} = diffusion - drift;
   
end


end