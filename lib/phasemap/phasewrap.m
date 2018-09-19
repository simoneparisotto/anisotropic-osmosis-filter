function Zw = phasewrap(Z,varargin) 
% phasewrap wraps values to the range -pi to pi or -180 to 180. 
%
%% Syntax
% 
%  Zw = phasewrap(Z) 
%  Zw = phasewrap(Z,'deg')  
%
%% Description 
% 
% Zw = phasewrap(Z) wraps angle Z in radians to the range -pi to pi. 
% 
% Zw = phasewrap(Z,'deg') wraps angle Z in degree to the range -180 to 180. 
% 
%% Example 
% 
% Z = peaks(900); 
% Zw = phasewrap(Z); 
% imagesc(Zw) 
% phasemap
% colorbar 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at Austin
% Institute for Geophysics (UTIG), May 2016. 
% http://www.chadagreene.com
% 
% If the phasemap function is useful for you, please do us a kindness and cite us: 
% 
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 
% 2016. True colors of oceanography: Guidelines for effective and accurate 
% colormap selection. Oceanography 29(3):9?13. 
% http://dx.doi.org/10.5670/oceanog.2016.66
% 
% See also phasemap, angle, and unwrap. 

%% Input checking: 

narginchk(1,2) 
assert(isnumeric(Z)==1,'Input Error: Z must be numeric.') 
if any(strncmpi(varargin,'deg',3))
   MaxAngle = 180; 
else
   MaxAngle = pi; 
end

%% Perform mathematics: 

Zw = mod(Z+MaxAngle,2*MaxAngle)-MaxAngle; 

end