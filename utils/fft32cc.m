function res = fft32cc(x)
% Use this function to fourier transform a bunch of 2D data sets at once

    xtemp = x(:,:,1);
    res   = 1/sqrt(length(xtemp(:)))*ifftshift(fft2(fftshift(fftshift(x,2)))); 

end

 