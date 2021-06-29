function [y_filt] = filter_kspace(y_raw,filt)
NX         = size(y_raw,1);
NY         = size(y_raw,2);
NZ         = size(y_raw,3);
axiX       = linspace(-NX/2,NX/2,NX);
axiY       = linspace(-NY/2,NY/2,NY);
axiZ       = linspace(-NZ/2,NZ/2,NZ);
filterVecX = 1-filt*(cos(0.5*pi*(axiX-NX/2)/(NX-NX/2))).^2;
filterVecY = 1-filt*(cos(0.5*pi*(axiY-NY/2)/(NY-NY/2))).^2;
filterVecZ = 1-filt*(cos(0.5*pi*(axiZ-NZ/2)/(NZ-NZ/2))).^2;
[x, y, z]  = meshgrid(filterVecY,filterVecX,filterVecZ);
filterMat  = x.*y.*z;
y_filt     = y_raw.*filterMat;
end

