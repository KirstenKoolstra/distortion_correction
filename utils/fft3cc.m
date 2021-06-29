function res = fft3cc(x)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

res = 1/sqrt(length(x(:)))*fftshift(fftn(ifftshift(x,1)));

end

  