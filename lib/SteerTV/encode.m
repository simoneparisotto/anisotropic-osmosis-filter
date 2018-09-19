%Author: Emmanuel Maggiori. March 2014.

%Input:
%origimage: input image

%Output:
%saliency and orientation fields of tensorized gradient

function [saliency,orientation] = encode(origimage) 


	[Dx,Dy] = gradient(double(origimage));

	saliency=double(abs(Dx)+abs(Dy));

    % take the orientation of the main direction
	orientation=pi/2+atan2(double(Dy),double(Dx));

end
