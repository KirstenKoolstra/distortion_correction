function res = ifft32cc(x)
% Use this function to fourier transform a bunch of 2D data sets at once

xtemp = x(:,:,1);
res   = sqrt(length(xtemp(:)))*ifftshift(ifftshift(ifft2(fftshift(x))),2);  

end

 