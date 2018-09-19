%Author: Emmanuel Maggiori. March 2014.

function kernel = w(m,h,w,sigma)

	w2=ceil(w/2.0);
	h2=ceil(h/2.0);



	[x, y] = meshgrid(1:w,1:h);

	kernel = exp(-(((x-w2).^2+(y-h2).^2)/(2*sigma^2))).*(((x-w2)+1i*(y-h2))./(sqrt((x-w2).^2+(y-h2).^2))).^m;

	kernel(h2,w2)=1;

end
