%% |phasebar| documentation
% The |phasebar| function places a circular donunt-shaped colorbar for phase 
% from -pi to pi or -180 degrees to 180 degrees. 
% 
%% Syntax
% 
%  phasebar
%  phasebar(...,'location',Location) 
%  phasebar(...,'size',Size) 
%  phasebar('deg') 
%  phasebar('rad') 
%  ax = phasebar(...) 
% 
%% Description 
% 
% |phasebar| places a donut-shaped colorbar on the current axes. 
%
% |phasebar(...,'location',Location)| specifies the corner (e.g., |'northeast'| or |'ne'|) 
% of the current axes in which to place the phasebar. Default location is the upper-right or |'ne'| 
% corner. 
%
% |phasebar(...,'size',Size)| specifies a size fraction of the current axes.  Default is 0.3. 
%
% |phasebar('deg')| plots labels at every 90 degrees. 
%
% |phasebar('rad')| plots labels at every pi/2 radians. 
%
% |ax = phasebar(...)| returns a handle |ax| of the axes in which the new axes are plotted. 
% 
%% Example

Z = 200*peaks(900); 
Zw = phasewrap(Z,'degrees'); 
imagesc(Zw) 
phasemap(12)
phasebar('location','se')

%% 
% Do you prefer a large phasebar with units of radians in the upper right-hand corner?  

phasebar('location','nw','size',0.5,'rad') 

%% 
% And instead of 12 discrete levels do you prefer a continuum? 

phasemap

%% Author Info
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the University of Texas 
% at Austin's Institute for Geophysics (UTIG), May 2016. 
% This function includes Kelly Kearney's |plotboxpos| function as a subfunction. 