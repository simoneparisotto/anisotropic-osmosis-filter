%% MATLAB Codes for the ANISOTROPIC OSMOSIS - REAL EXAMPLE
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

%% PARAMETERS
% time domain after which I update W and go on till tt(end)
time         = [0 1000];
MAXITER      = 100;
printstep    = 1000;

imagenumber = [6 14 15 44];

%% DO NOT TOUCH FROM HERE
reducecolorimage = false;

compute_osmosis  = true;
show_osmosis     = false;

flag_classic    = 1;
flag_classic_w  = 0;
flag_coherence  = 0;
flag_mirebeau   = 1;


%% MAIN LOOP
% parpool(min(maxNumCompThreads,5));
%parpool(6);
%parfor imn = 1:numel(imagenumber)
for imn = 1:numel(imagenumber)
    
    caseimage = imagenumber(imn);
    
    % CREATE DIRECTORIES
    dir_result = ['./results/',num2str(caseimage)];
    dir_result_mat = ['./results/mat/'];
    
    if ~exist(dir_result,'dir')
        mkdir(dir_result);
    end
    if ~exist(dir_result_mat,'dir')
        mkdir(dir_result_mat);
    end
    
    % READ IMAGE AND MASK + sigma,rho for structure tensor
    imagefile = ['./../dataset/shadow/',num2str(caseimage),'.png'];
    maskfile  = ['./../dataset/mask/',num2str(caseimage),'.png'];
    
    umat     = im2double(imread(imagefile));
    umask    = repmat(sum(im2double(imread(maskfile)),3)>0,[1 1 size(umat,3)]);
    
    if reducecolorimage && size(umat,3)>1
        umat = rgb2gray(umat);
    end
    
    pos_offset    = 1;
    v             = umat+pos_offset; % is the image for d=\grad\log v
    if flag_classic
        u_CLASSIC     = umat+pos_offset; % image to evolve in time
    end
    if flag_classic_w
        u_CLASSIC_W   = umat+pos_offset; % image to evolve in time
    end
    if flag_coherence
        u_COHERENCE   = umat+pos_offset; % image to evolve in time
    end
    if flag_mirebeau
        u_MIREBEAU    = umat+pos_offset; % image to evolve in time
    end
    
    %% SPATIAL DISCRETIZATION
    P3classic   = cell(3,1);
    P3classic_W = cell(3,1);
    P3coherence = cell(3,1);
    P3mirebeau  = cell(3,1);
    
    local_load = load(['./results/mat/',num2str(caseimage),'.mat'],'TVF','b','V','W');
    V  = local_load.V.*(1-repmat(umask(:,:,1),1,1,2));
    W  = local_load.W;
    A  = local_load.TVF.A;
    b  = local_load.b;
    
    % SPACE DISCRETISATION
    for kk=1:size(umat,3)
        if flag_classic
            [~, ~, P3classic{kk}] = classical_osmosis_discretization(v(:,:,kk),umask(:,:,kk));
        end
        if flag_classic_w
            P3classic_W{kk}       = classical_osmosis_discretization_W(v(:,:,kk),umask(:,:,kk),W);
        end
        if flag_coherence
            P3coherence{kk}       = coherence_osmosis_discretization(v(:,:,kk),umask(:,:,kk),W);
        end
        if flag_mirebeau
            flag_upwind  = 0;
            P3mirebeau{kk}        = create_mirebeau_matrix(v(:,:,kk),umask(:,:,kk),W,flag_upwind);
        end
    end
    
    %% TIME INTEGRATION WITH DIFFERENT SHADOW REMOVAL METHODS
    for tt = 1:MAXITER
        
        for kk=1:size(umat,3)
            if flag_classic
                u_CLASSIC(:,:,kk)   = reshape(expleja(time(2),P3classic{kk}{1},reshape(u_CLASSIC(:,:,kk),[],1)),size(u_CLASSIC,1),size(u_CLASSIC,2));
            end
            if flag_classic_w
                u_CLASSIC_W(:,:,kk) = reshape(expleja(time(2),P3classic_W{kk}{1},reshape(u_CLASSIC_W(:,:,kk),[],1)),size(u_CLASSIC_W,1),size(u_CLASSIC_W,2));
            end
            if flag_coherence
                u_COHERENCE(:,:,kk) = reshape(expleja(time(2),P3coherence{kk}{1},reshape(u_COHERENCE(:,:,kk),[],1)),size(u_COHERENCE,1),size(u_COHERENCE,2));
            end
            if flag_mirebeau
                u_MIREBEAU(:,:,kk)  = reshape(expleja(time(2),P3mirebeau{kk},reshape(u_MIREBEAU(:,:,kk),[],1)),size(u_MIREBEAU,1),size(u_MIREBEAU,2));
            end
        end
        
        
        %% FIGURE
        if show_osmosis
            
            [X,Y] = meshgrid(1:size(umat,2),size(umat,1):-1:1);
            
            h1=figure(1);
            
            subplot(2,4,1)
            imshow(v-pos_offset,[0,1])
            title(['T = 0'])
            
            
            subplot(2,4,2)
            imagesc(X(:),Y(:),v-pos_offset)
            if size(v,3)==1
                colormap(gray(256))
            end
            axis off
            axis image
            hold on
            quiver(X,Y,V(:,:,1),V(:,:,2));
            title('vector field z perp')
            hold off
            set(gca,'YDir','normal')
            
            
            subplot(2,4,3)
            imshow(b{1}.*(1-umask(:,:,1)),[0,1])
            title('Anisotropy b1')
            
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
            title(['COHERENCE, T = ',num2str(tt*time(2))])
            end
            if flag_mirebeau
            subplot(2,4,8)
            imshow(u_MIREBEAU-pos_offset,[0,1])
            title(['MIREBEAU, T = ',num2str(tt*time(2))])
            end
            
            truesize(h1,[600 600])
            
            pause(0.01)
            
        end
        
        % SAVE FIGURES at timestep tt*time(2)
%         if ~mod(tt*time(2),time(2))
%             if flag_classic
%                 imwrite(u_CLASSIC-pos_offset,  [dir_result,'/',num2str(caseimage),'_classic_',num2str(tt*time(2)),'.png'])
%             end
%             if flag_classic_w
%                 imwrite(u_CLASSIC_W-pos_offset,[dir_result,'/',num2str(caseimage),'_classic_W_',num2str(tt*time(2)),'.png'])
%             end
%             if flag_coherence
%                 imwrite(u_COHERENCE-pos_offset,[dir_result,'/',num2str(caseimage),'_coherence_',num2str(tt*time(2)),'.png'])
%             end
%             if flag_mirebeau
%                 imwrite(u_MIREBEAU-pos_offset, [dir_result,'/',num2str(caseimage),'_mirebeau_',num2str(tt*time(2)),'.png'])
%             end
%         end
        
    end
    
    if flag_classic
        imwrite(u_CLASSIC-pos_offset,  [dir_result,'/',num2str(caseimage),'_classic.png'])
    end
    if flag_classic_w
        imwrite(u_CLASSIC_W-pos_offset,[dir_result,'/',num2str(caseimage),'_classic_W.png'])
    end
    if flag_coherence
        imwrite(u_COHERENCE-pos_offset,[dir_result,'/',num2str(caseimage),'_coherence.png'])
    end
    if flag_mirebeau
        imwrite(u_MIREBEAU-pos_offset, [dir_result,'/',num2str(caseimage),'_mirebeau.png'])
    end
    
end


%delete(gcp('nocreate'))