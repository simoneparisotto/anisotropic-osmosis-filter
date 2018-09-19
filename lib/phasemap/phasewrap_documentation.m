%% |phasewrap| documentation 
% The |phasewrap| function wraps values to the range -pi to pi or -180 to 180. 
%
%% Syntax
% 
%  Zw = phasewrap(Z) 
%  Zw = phasewrap(Z,'deg')  
%
%% Description 
% 
% |Zw = phasewrap(Z)| wraps angle |Z| in radians to the range -pi to pi. 
% 
% |Zw = phasewrap(Z,'deg')| wraps angle |Z| in degree to the range -180 to 180. 
% 
%% Example 

Z = 2*peaks(900); 
Zw = phasewrap(Z); 
imagesc(Zw) 
phasemap
cb = colorbar; 
ylabel(cb,' phase ') 

%% Author Info
% This function was written by Chad A. Greene of the University of Texas at Austin
% Institute for Geophysics (UTIG), May 2016. 
% <http://www.chadagreene.com>