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
%  April, 2018
%%

clear
close all
clc

addpath('./../dataset/')
addpath('./../lib/')
addpath('./../lib/expleja/')
addpath('./../lib/Mirebeau')
addpath('./../lib/Mirebeau/ParallelMatrix')
addpath('./../lib/phasemap')
addpath('./../lib/SteerTV/')
addpath('./../lib/export_fig-master')

%% PARAMETERS
% b is the minor axis of ellipse b=1 is the circle b=0 is the segment (degenerate ellipse)
b{1}  = 0.05;
b{2}  = 1;
sigma = 0.25;
rho   = 2;

imagenumber = [5 6 14 15 44 66];

%% DO NOT TOUCH FROM HERE
reducecolorimage = false;

redoexperiment   = true;
compute_fields   = true;
print_fields     = true;


%% MAIN LOOP

for imn = 1:numel(imagenumber)
    
    % CREATE DIRECTORIES
    dir_result     = ['./results/',num2str(imagenumber(imn))];
    dir_result_mat = ['./results/mat/'];
    
    if ~exist(dir_result,'dir')
        mkdir(dir_result);
    end
    if ~exist(dir_result_mat,'dir')
        mkdir(dir_result_mat);
    end
    
    % READ IMAGE AND MASK + sigma,rho for structure tensor
    imagefile = ['./../dataset/shadow/',num2str(imagenumber(imn)),'.png'];
    maskfile  = ['./../dataset/mask/',num2str(imagenumber(imn)),'.png'];
    
    umat     = im2double(imread(imagefile));
    umask    = repmat(sum(im2double(imread(maskfile)),3)>0,[1 1 size(umat,3)]);
   
    if reducecolorimage && size(umat,3)>1
        umat = rgb2gray(umat);
    end
    
    %% CREATE V AND W
    if  ~exist([dir_result_mat,num2str(imagenumber(imn)),'.mat'],'file') || redoexperiment
        [STF,TVF] = computefield(umat,umask(:,:,1),sigma,rho);
        save([dir_result_mat,num2str(imagenumber(imn)),'.mat'],'STF','TVF','b');
    else
        load([dir_result_mat,num2str(imagenumber(imn)),'.mat'],'STF','TVF','b');
    end
    
    if compute_fields
        % rotate z by 90 to match the discretization used
        [V,W] = compute_W(TVF.z,b,umask);
        save([dir_result_mat,num2str(imagenumber(imn)),'.mat'],'V','W','-append')
    end
    
    %% PLOT
    if print_fields
        
        % This is the cartesian grid: (1,1) bottom-left, (M,N) top-right
        [X,Y] = meshgrid(1:size(umat,2),size(umat,1):-1:1);
        
        step    = 10;
        scaleph = 0.25;
        if imagenumber(imn)==14
            posph = 'ne';
        else
            posph   = 'sw';
        end
       
        % rotate to match the image axes
        TVF_z_plot      = TVF.z_perp;
        TVF_z_perp_plot = TVF.z;
        STF_z_plot      = STF.z_perp;
        STF_z_perp_plot = STF.z;
        
        % angles
        theta_z_plot         = TVF.theta_z;
        theta_z_perp_plot    = TVF.theta_z_perp;
        theta_z_mask         = TVF.theta_z;
        theta_z_perp_mask    = TVF.theta_z_perp;
        theta_z_mask(umask(:,:,1)>0)      = NaN;
        theta_z_perp_mask(umask(:,:,1)>0) = NaN;
        
        theta_zold_plot      = STF.theta_z;
        theta_zold_perp_plot = STF.theta_z_perp;
        theta_zold_mask      = STF.theta_z;
        theta_zold_perp_mask = STF.theta_z_perp;
        theta_zold_mask(umask(:,:,1)>0)      = NaN;
        theta_zold_perp_mask(umask(:,:,1)>0) = NaN;
        
        %%
        h_z = figure(101);
        imagesc(X(:),Y(:),umat),
        if size(umat,3)==1
            colormap(gray(256))
        end
        axis off
        axis image
        hold on
        streamline(X,Y,TVF_z_plot(:,:,1),TVF_z_plot(:,:,2),X(1:step:end,1:step:end),Y(1:step:end,1:step:end))
        hold off
        set(gca,'YDir','normal')
        
        h_z_perp = figure(102);
        imagesc(X(:),Y(:),umat)
        if size(umat,3)==1
            colormap('gray')
        end
        axis off
        axis image
        hold on
        streamline(X,Y,TVF_z_perp_plot(:,:,1),TVF_z_perp_plot(:,:,2),X(1:step:end,1:step:end),Y(1:step:end,1:step:end))
        hold off
        set(gca,'YDir','normal')
        
        h_zold = figure(103);
        imagesc(X(:),Y(:),umat)
        if size(umat,3)==1
            colormap(gray(256))
        end
        axis off
        axis image
        hold on
        streamline(X,Y,STF_z_plot(:,:,1),STF_z_plot(:,:,2),X(1:step:end,1:step:end),Y(1:step:end,1:step:end))
        hold off
        set(gca,'YDir','normal')
        
        h_zold_perp = figure(104);
        imagesc(X(:),Y(:),umat)
        if size(umat,3)==1
            colormap(gray(256))
        end
        axis off
        axis image
        hold on
        streamline(X,Y,STF_z_perp_plot(:,:,1),STF_z_perp_plot(:,:,2),X(1:step:end,1:step:end),Y(1:step:end,1:step:end))
        hold off
        set(gca,'YDir','normal')
        
        h_l1 = figure(105);
        imagesc(X(:),Y(:),TVF.l1)
        colormap('gray')
        axis off
        axis image
        set(gca,'YDir','normal')
        
        h_l1old = figure(106);
        imagesc(X(:),Y(:),STF.l1)
        colormap('gray')
        axis off
        axis image
        set(gca,'YDir','normal')
        
        % PLOT ANGLES
        h_theta_z = figure(207);
        imagesc(X(:),Y(:),theta_z_plot);
        axis off
        axis image
        phasemap(12);
        phasebar('location',posph,'size',scaleph)
        set(gca,'YDir','normal')
        
        h_theta_z_perp = figure(107);
        imagesc(X(:),Y(:),theta_z_perp_plot);
        axis off
        axis image
        phasemap(12);
        phasebar('location',posph,'size',scaleph) 
        set(gca,'YDir','normal')
                
        h_theta_zold = figure(208);
        imagesc(X(:),Y(:),theta_zold_plot);
        axis off
        axis image
        phasemap(12);
        phasebar('location',posph,'size',scaleph)
        set(gca,'YDir','normal')
        
        h_theta_zold_perp = figure(108);
        imagesc(X(:),Y(:),theta_zold_perp_plot);
        axis off
        axis image
        phasemap(12);
        phasebar('location',posph,'size',scaleph)
        set(gca,'YDir','normal')
             
        h_thetamask = figure(109);
        imagesc(X(:),Y(:),theta_z_perp_mask,'AlphaData',1-umask(:,:,1));
        axis image
        set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor',[1,1,1],'YColor',[1,1,1],'color',[1 1 1]);
        phasemap(180);
        phasebar('location',posph,'size',scaleph)
        set(gca,'YDir','normal')
        
        h_thetaoldmask = figure(110);
        imagesc(X(:),Y(:),theta_zold_perp_mask,'AlphaData',1-umask(:,:,1));
        axis image
        set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'XColor',[1,1,1],'YColor',[1,1,1],'color',[1 1 1]);
        phasemap(180);
        phasebar('location',posph,'size',scaleph)
        set(gca,'YDir','normal')
        
        
        %% SAVE FIGURES
        export_fig(h_z,              [dir_result,'/',num2str(imagenumber(imn)),'_TVF_z','.png'], '-m3')
        export_fig(h_z_perp,         [dir_result,'/',num2str(imagenumber(imn)),'_TVF_z_perp','.png'], '-m3')
        export_fig(h_l1,             [dir_result,'/',num2str(imagenumber(imn)),'_TVF_l1','.png'], '-m3')
        export_fig(h_theta_z,        [dir_result,'/',num2str(imagenumber(imn)),'_TVF_z_theta','.png'], '-m3')
        export_fig(h_theta_z_perp,   [dir_result,'/',num2str(imagenumber(imn)),'_TVF_z_perp_theta','.png'], '-m3')
        
        export_fig(h_zold,           [dir_result,'/',num2str(imagenumber(imn)),'_STF_z','.png'], '-m3')
        export_fig(h_zold_perp,      [dir_result,'/',num2str(imagenumber(imn)),'_STF_z_perp','.png'], '-m3')
        export_fig(h_l1old,          [dir_result,'/',num2str(imagenumber(imn)),'_STF_l1','.png'], '-m3')
        export_fig(h_theta_zold,     [dir_result,'/',num2str(imagenumber(imn)),'_STF_z_theta','.png'], '-m3')
        export_fig(h_theta_zold_perp,[dir_result,'/',num2str(imagenumber(imn)),'_STF_z_perp_theta','.png'], '-m3')
        
        export_fig(h_thetamask,      [dir_result,'/',num2str(imagenumber(imn)),'_TVFmask_z_perp','.png'], '-m3')
        export_fig(h_thetaoldmask,   [dir_result,'/',num2str(imagenumber(imn)),'_STFmask_z_perp','.png'], '-m3')
        
        imwrite(b{1},[dir_result,'/',num2str(imagenumber(imn)),'_b1','.png'])
        imwrite(mat2gray(STF.A,[0,1]),[dir_result,'/',num2str(imagenumber(imn)),'_STF_A','.png'])
        imwrite(mat2gray(TVF.A,[0,1]),[dir_result,'/',num2str(imagenumber(imn)),'_TVF_A','.png'])
        
        close all
        
    end
    
end
