function [u,ushadow,umask,V,WW,WA,V_all] = create_synthetic(theta,b,select)

% '1' rotated greyscale stripes
% '2' bended greyscale stripes
% '3' bended rgb stripes
% '4' bamboo greyscale image, fixed 45 direction

% v is the tangential direction

coeff_shadow = 0.6;

switch select
    
    case '1' 
        stripes = repmat(cat(3,...
            [100*ones(101,20),   1*ones(101,10),200*ones(101,15),150*ones(101,5), 90*ones(101,25),   5*ones(101,10), 200*ones(101,5), 100*ones(101,11)],...
            [250*ones(101,20), 200*ones(101,10),  1*ones(101,15), 20*ones(101,5), 10*ones(101,25), 100*ones(101,10),  75*ones(101,5), 220*ones(101,11)],...
            [ 80*ones(101,20), 200*ones(101,10),  1*ones(101,15), 10*ones(101,5), 90*ones(101,25), 150*ones(101,10), 150*ones(101,5),  90*ones(101,11)]),[5,5,1]);
        stripes = stripes(:,:,1).'/255;
        
        u = rotateAround(stripes,253,253,theta,'bilinear');
        
        %u = imrotate(stripes,theta,'bilinear');
        u = u(253-100:253+100,253-100:253+100,:);
        
        % get angle of main direction in degree in (x,y) coordinates
        thetarad = deg2rad(theta);
        
        % get gradient in (i,j) coordinates
        V  = cat(3,-sin(thetarad+pi/2)*ones(size(u)),cos(thetarad+pi/2)*ones(size(u)));

        mask                  = ones(size(u));
        mask(90:160,90:160,:) = coeff_shadow;
        mask = imgaussfilt(mask,2);
        
        umask                = ones(size(mask));
        umask(85:165,85:165) = 0;
        umask(95:155,95:155) = 1;
        
        ushadow = mask .* u;
        
        
        
     case '11' 
         
        stripes = repmat(cat(3,...
            [100*ones(101,20),   1*ones(101,10),200*ones(101,15),150*ones(101,5), 90*ones(101,25),   5*ones(101,10), 200*ones(101,5), 100*ones(101,11)],...
            [250*ones(101,20), 200*ones(101,10),  1*ones(101,15), 20*ones(101,5), 10*ones(101,25), 100*ones(101,10),  75*ones(101,5), 220*ones(101,11)],...
            [ 80*ones(101,20), 200*ones(101,10),  1*ones(101,15), 10*ones(101,5), 90*ones(101,25), 150*ones(101,10), 150*ones(101,5),  90*ones(101,11)]),[5,5,1]);
        stripes = stripes(:,:,1).'/255;
        
        u = rotateAround(stripes,253,253,theta,'bilinear');
        
        %u = imrotate(stripes,theta,'bilinear');
        u = u(253-100:253+100,253-100:253+100,:);
        
        % get angle of main direction in degree in (x,y) coordinates
        thetarad = deg2rad(theta);
        
        % get gradient in (i,j) coordinates
        V  = cat(3,-sin(thetarad+pi/2)*ones(size(u)),cos(thetarad+pi/2)*ones(size(u)));

        mask = double(rot90(double(imread('../dataset/mask/hand.png'))))/255;
        mask(mask==0) = coeff_shadow;
        mask = imresize(mask,size(u));
        
        umask = imerode(1-edge(mask),ones(7));
        mask = imgaussfilt(mask,1.3);
        
        ushadow = mask .* u;
        
        
    case '2' 
        
        [X,Y] = ndgrid((1:201)/201,(1:201)/201);
        u = exp(-X.^2-Y.^2); %u = imnoise(u(:,:,1),'gaussian');
        u = mat2gray(u);
        V = grad_centered(u);
        V = cat(3,-V(:,:,2),V(:,:,1)); 
        normV = repmat(sqrt(sum(V.^2,3)),[1 1 2]);
        V = V./normV;

        mask                  = ones(size(u));
        mask(30:100,30:100,:) = coeff_shadow;
        mask = imgaussfilt(mask,2);
        
        umask                 = ones(size(mask));
        umask(25:105,25:105)  = 0;
        umask(35:95,35:95)  = 1;
        
        ushadow = round(5*u)/5;
        
        ushadow = mask .* ushadow;
        
    case '3'
        
        [X,Y] = ndgrid((1:201)/201,(1:201)/201);
        u = exp(-X.^2-Y.^2); %u = imnoise(u(:,:,1),'gaussian');
        u = mat2gray(u);
        V = grad_centered(u);
        %V = cat(3,-V(:,:,2),V(:,:,1)); 
        normV = repmat(sqrt(sum(V.^2,3)),[1 1 2]);
        V = V./normV;


        mask                  = ones(size(u));
        mask(30:100,30:100,:) = coeff_shadow;
        mask = imgaussfilt(mask,1.3);
        
        umask                 = ones(size(mask));
        umask(25:105,25:105)  = 0;
        umask(35:95,35:95)  = 1;
         
        ushadow = round(5*u)/5;         
        
        nlevels = unique(ushadow(:));
        colors  = jet(numel(nlevels));
        masklev = zeros(size(ushadow));
        for kk = 1:numel(nlevels)
            ulev{kk} = cat(3,colors(kk,1)*ones(size(ushadow)),colors(kk,2)*ones(size(ushadow)),colors(kk,3)*ones(size(ushadow)));
            masklev(ushadow==nlevels(kk)) = kk;
        end
        masklev = repmat(masklev,1,1,3);
        ucolor = zeros([size(ushadow),3]);
        for kk=1:numel(nlevels)
            ucolor(masklev==kk) = ulev{kk}(masklev==kk);
        end
        
        mask = repmat(mask,1,1,3);
        
        %ucolor(mask==0) = 0.7*ucolor(mask==0);
        ucolor = mask .* ucolor;
        ushadow = ucolor;
        
    case '33'
        
        [X,Y] = ndgrid((1:201)/201,(1:201)/201);
        u = exp(-X.^2-Y.^2); %u = imnoise(u(:,:,1),'gaussian');
        u = mat2gray(u);
        V = grad_centered(u);
        %V = cat(3,-V(:,:,2),V(:,:,1)); 
        normV = repmat(sqrt(sum(V.^2,3)),[1 1 2]);
        V = V./normV;


        mask = double(flipud(rot90(rot90(rot90(double(imread('../dataset/mask/hand.png')))))))/255;
        mask(mask==0) = coeff_shadow;
        mask = imresize(mask,size(u));
        
        umask = imerode(1-edge(mask),ones(7));
        mask  = imgaussfilt(mask,1.3);
         
        ushadow = round(5*u)/5;         
        
        nlevels = unique(ushadow(:));
        colors  = jet(numel(nlevels));
        masklev = zeros(size(ushadow));
        for kk = 1:numel(nlevels)
            ulev{kk} = cat(3,colors(kk,1)*ones(size(ushadow)),colors(kk,2)*ones(size(ushadow)),colors(kk,3)*ones(size(ushadow)));
            masklev(ushadow==nlevels(kk)) = kk;
        end
        masklev = repmat(masklev,1,1,3);
        ucolor = zeros([size(ushadow),3]);
        for kk=1:numel(nlevels)
            ucolor(masklev==kk) = ulev{kk}(masklev==kk);
        end
        
        mask = repmat(mask,1,1,3);
        
        %ucolor(mask==0) = 0.7*ucolor(mask==0);
        ucolor = mask .* ucolor;
        ushadow = ucolor;
        
    case '4'
        u        = im2double(imread('./../dataset/shadow/4.png'));
        %ushadow  = im2double(imread('./../dataset/shadow/4.png'));
        %umask    = im2double(imread('./../dataset/mask/4.png'));
        
        
        mask                  = ones(size(u));
        mask(100:200,100:200,:) = coeff_shadow;
        mask = imgaussfilt(mask,2);
        
        umask                   = ones(size(mask));
        umask(95:205,95:205)  = 0;
        umask(105:195,105:195)  = 1;
        
        ushadow = mask .* u;
        
        Vij = cat(3,...
            -sin(deg2rad(45)).*ones(size(u,1),size(u,2)),...
            cos(deg2rad(45)).*ones(size(u,1),size(u,2)));
        
        Vxy = cat(3,...
            cos(deg2rad(45)).*ones(size(u,1),size(u,2)),...
            sin(deg2rad(45)).*ones(size(u,1),size(u,2)));
        
        V = Vij;


        
end

B1 = b(1);
B2 = b(2);

W1 = cat(3,ones(size(u)),zeros(size(u)),zeros(size(u)),ones(size(u)));
WW = cat(3, B1.^2.*V(:,:,1).^2         + B2.^2.*V(:,:,2).^2,...
            B1.^2.*V(:,:,1).*V(:,:,2)  - B2.^2.*V(:,:,1).*V(:,:,2),...
            B1.^2.*V(:,:,1).*V(:,:,2)  - B2.^2.*V(:,:,1).*V(:,:,2),...
            B1.^2.*V(:,:,2).^2         + B2.^2.*V(:,:,1).^2);
WA = WW; % save before modifying (just in case)
WW(repmat(umask,1,1,4)==1) = W1(repmat(umask,1,1,4)==1);
V_all = V;
V(repmat(umask,1,1,2)==1)  = 0;