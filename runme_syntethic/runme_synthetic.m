%% MATLAB Codes for the ANISOTROPIC OSMOSIS - SYNTHETIC EXAMPLE
%  Copyright (c) 2018, Simone Parisotto
%  All rights reserved.
%
%  Author:
%  Simone Parisotto (email: sp751 at cam dot ac dot uk)
%
%  Address:
%  Cambridge Image Analysis
%  Centre for Mathematical Sciences
%  Wilberforce Road
%  Cambridge CB3 0WA
%  United Kingdom
%
%  Date:
%  September, 2018
%%

clear
close all
clc

addpath('./../dataset/')
addpath('./../lib/')
addpath('./../lib/operators')
addpath('./../lib/Mirebeau')
addpath('./../lib/Mirebeau/ParallelMatrix')
addpath('./../lib/SteerTV/')
addpath('./../lib/expleja/')
addpath('./../lib/phasemap')
addpath('./../lib/export_fig-master')

%% IN OUR CODE
% * b1 is the minor axis of ellipse b1=1 is the circle b1=0 is the
%  segment  (degenerate ellipse)
% * v has been choosen as the gradient direction
% * b = [b1 b2] where b1 is for v and b2 is for v orthogonal
% * umask = 1 then: b=[1,1];
% * umask = 0 then: b=[b1,1];
b = [0.05,1];

%% LOAD CASE:
% '1' rotated greyscale stripes
% '2' bended greyscale stripes
% '3' bended rgb stripes
% '4' bamboo greyscale image, fixed 45 direction
% '11' rotated greyscale stripes with hand shandow
% '33' bended rgb stripes with hand shandow
imagenumber = [1,3,4,11,33];
imagenumber = [11,33];

%% SELECT SPACE DISCRETISATION
flag_classic    = 1; % CLASSIC OSMOSIS FILTER
flag_classic_w  = 0; % ALTERNATIVE OF CLASSIC WITH ANISOTROPIC WEIGHTS
flag_coherence  = 0; % COHERENCE ENHANCING
flag_mirebeau   = 1; % SPARSE NON-NEGATIVE STENCIL (MIREBEAU)


%% DO NOT TOUCH FROM HERE

step_plot = 1000;

for imn = 1:numel(imagenumber)
    
    caseimage = num2str(imagenumber(imn));
    switch caseimage
        case {'1','11'}
            step = 5;
            angle = (0:step:360-step);
            angle = 65;
        otherwise
            angle = NaN;
    end
    
    if ~exist(['results/case',caseimage],'dir')
        mkdir(['results/case',caseimage])
    end
    pos_offset = 1;
    
    for aa = 1:numel(angle) % IN CASE OF MYULTIPLE DIRECTION TESTS
        
        dir_result = ['./results/case',caseimage];
        
        [ureal,umat,umask,V,W,WA,V_all] = create_synthetic(angle(aa),b,caseimage);
        
        v             = umat+pos_offset; % is the image for d=\grad\log v
        
        if flag_classic
            u_CLASSIC     = umat+pos_offset; % image to evolve in time
        end
        
        if flag_classic_w
            u_CLASSIC_W   = umat+pos_offset; % image to evolve in time
        end
        
        if flag_coherence
            u_COHERENCE = umat+pos_offset; % image to evolve in time
        end
        
        if flag_mirebeau
            u_MIREBEAU    = umat+pos_offset; % image to evolve in time
        end
        
        % time domain after which I update W and go on till tt(end)
        time     = [0 100];
        timespan = 1:100;
        
        %% PLOT
        
        [X,Y] = meshgrid(1:size(umat,2),size(umat,1):-1:1);
        step    = 10;
        scaleph = 0.25;
        posph   = 'sw';
        
        h_e1 = figure(101);
        imagesc(X(:),Y(:),umat)
        if size(umat,3)==1
            colormap(gray(256))
        end
        axis off
        hold on
        axis image
        % gradient in (x,y) since V_all is in (i,j)
        streamline(X,Y,V_all(:,:,2),-V_all(:,:,1),X(1:step:end,1:step:end),Y(1:step:end,1:step:end))
        hold off
        set(gca,'YDir','normal')
        
        h_e2 = figure(102);
        imagesc(X(:),Y(:),umat)
        if size(umat,3)==1
            colormap(gray(256))
        end
        axis off
        hold on
        axis image
        % directional in (x,y) since V_all is in (i,j)
        streamline(X,Y,V_all(:,:,1),V_all(:,:,2),X(1:step:end,1:step:end),Y(1:step:end,1:step:end))
        hold off
        set(gca,'YDir','normal')
        
        % create angle of main direction in (x,y) coordinates
        % remember that V is the gradient in (i,j) so we need
        % firstly to project on (x,y) and then to take the orthogonal
        theta     = (atan2d(-V_all(:,:,1),V_all(:,:,2)))-90;
        thetamask = theta;
        thetamask(umask(:,:,1)>0) = NaN;
        
        h_thetamask = figure(109);
        imagesc(X(:),Y(:),thetamask,'AlphaData',1-umask(:,:,1));
        set(gca,'YDir','normal')
        axis image
        set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor',[1,1,1],'YColor',[1,1,1],'color',[1 1 1]);
        phasemap(12,'degrees');
        phasebar('location',posph,'size',scaleph)
        
        export_fig(h_e1,[dir_result,'/v',num2str(angle(aa)),'.png'], '-m3')
        export_fig(h_e2,[dir_result,'/v',num2str(angle(aa)),'_perp.png'], '-m3')
        export_fig(h_thetamask,[dir_result,'/theta',num2str(angle(aa)),'.png'], '-m3')
        
        
        %% ANISOTROPIC SHADOW REMOVAL
        
        for kk=1:size(umat,3)
            
            %% OPERATORS FOR SPATIAL DISCRETIZATION
            % compute the matrices only once, since we exactly know the direction
            if flag_classic
                [~, ~, P3classic{kk}] = classical_osmosis_discretization(v(:,:,kk),umask);
            end
            if flag_classic_w
                P3classic_W{kk}       = classical_osmosis_discretization_W(v(:,:,kk),umask,W);
            end
            if flag_coherence
                P3coherence{kk}       = coherence_osmosis_discretization(v(:,:,kk),umask,W);
            end
            if flag_mirebeau
                flag_upwind = 1;
                P3mirebeau{kk}        = create_mirebeau_matrix(v(:,:,kk),umask,W,flag_upwind);
            end
        end
        
        %% SOLVE IN TIME
        for tt = timespan
            for kk=1:size(umat,3)
                
                %% CLASSIC SHADOW REMOVAL
                if flag_classic
                    u_CLASSIC(:,:,kk) = reshape(expleja(time(2),P3classic{kk}{1},reshape(u_CLASSIC(:,:,kk),[],1)),size(u_CLASSIC,1),size(u_CLASSIC,2));
                end
                
                %% CLASSIC WITH W SHADOW REMOVAL
                if flag_classic_w
                    u_CLASSIC_W(:,:,kk) = reshape(expleja(time(2),P3classic_W{kk}{1},reshape(u_CLASSIC_W(:,:,kk),[],1)),size(u_CLASSIC_W,1),size(u_CLASSIC_W,2));
                end
                
                %% COHERENCE SCHEME SHADOW REMOVAL
                if flag_coherence
                    u_COHERENCE(:,:,kk) = reshape(expleja(time(2),P3coherence{kk}{1},reshape(u_COHERENCE(:,:,kk),[],1)),size(u_COHERENCE,1),size(u_COHERENCE,2));
                end
                
                %% MIREBEAU SHADOW REMOVAL
                if flag_mirebeau
                    u_MIREBEAU(:,:,kk) = reshape(expleja(time(2),P3mirebeau{kk},reshape(u_MIREBEAU(:,:,kk),[],1)),size(u_MIREBEAU,1),size(u_MIREBEAU,2));
                end
                
            end
            
            % PLOT RESULTS
            if ~mod(time(2),step_plot)
                figure(1)
                subplot(2,4,1)
                imshow(v-pos_offset,[0,1])
                title(['T = 0'])
                subplot(2,4,2)
                imshow(v-pos_offset,[0,1])
                hold on
                quiver(V(:,:,2),V(:,:,1));
                title(['T = ',num2str(tt),' vector field'])
                hold off
                
                if flag_classic
                    subplot(2,4,5)
                    imshow(u_CLASSIC-pos_offset,[0,1])
                    title(['CLASSIC, T = ',num2str(tt*time(2))])
                end
                
                if flag_classic_w
                    subplot(2,4,6)
                    imshow(u_CLASSIC_W-pos_offset,[0,1])
                    title(['CLASSIC with W, T = ',num2str(tt*time(2))])
                end
                
                if flag_coherence
                    subplot(2,4,7)
                    imshow(u_COHERENCE-pos_offset,[0,1])
                    title(['COHERENCE, b=(',num2str(b(1)),',',num2str(b(2)),'), T = ',num2str(tt*time(2))])
                end
                
                if flag_mirebeau
                    subplot(2,4,8)
                    imshow(u_MIREABEAU-pos_offset,[0,1])
                    title(['MIREBEAU, b=(',num2str(b(1)),',',num2str(b(2)),'), T = ',num2str(tt*time(2))])
                end
                pause(0.01)
            end
            
        end
        
        % WRITE RESULTS
        switch caseimage
            
            case {'1','11'}
                imwrite(umat,[dir_result,'/u_start_',num2str(angle(aa)),'.png']);
                imwrite(umask,[dir_result,'/u_mask_',num2str(angle(aa)),'.png']);
                if flag_classic
                    imwrite(u_CLASSIC-pos_offset,[dir_result,'/u_CLA_',num2str(angle(aa)),'.png']);
                end
                if flag_classic_w
                    imwrite(u_CLASSIC_W-pos_offset,[dir_result,'/u_CLA_W_',num2str(angle(aa)),'.png']);
                end
                if flag_coherence
                    imwrite(u_COHERENCE-pos_offset,[dir_result,'/u_DIR_',num2str(angle(aa)),'.png']);
                end
                if flag_mirebeau
                    imwrite(u_MIREBEAU-pos_offset,[dir_result,'/u_MIR_',num2str(angle(aa)),'.png']);
                end
                
            case {'2','3','33','4'}
                imwrite(umat,[dir_result,'/u_start.png']);
                imwrite(umask,[dir_result,'/u_mask.png']);
                if flag_classic
                    imwrite(u_CLASSIC-pos_offset,[dir_result,'/u_CLA.png']);
                end
                if flag_classic_w
                    imwrite(u_CLASSIC_W-pos_offset,[dir_result,'/u_CLA_W','.png']);
                end
                if flag_coherence
                    imwrite(u_COHERENCE-pos_offset,[dir_result,'/u_DIR.png']);
                end
                if flag_mirebeau
                    imwrite(u_MIREBEAU-pos_offset,[dir_result,'/u_MIR.png']);
                end
                
        end
        
    end
end