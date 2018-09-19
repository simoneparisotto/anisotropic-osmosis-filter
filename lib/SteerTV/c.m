%Author: Emmanuel Maggiori. March 2014.

function im = c(m,s,be) %s: stickness field, be:orientation field

	im=s.*exp(-1i*m*be);


end
