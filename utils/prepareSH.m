function [spherHarm_res] = prepareSH(nx,ny)

fieldMapResolution = 5;

coordDimX = linspace(-fieldMapResolution*(nx-1)/2,fieldMapResolution*(nx-1)/2, nx);
coordDimY = linspace(-fieldMapResolution*(ny-1)/2,fieldMapResolution*(ny-1)/2, ny);

[coordx, coordy] = meshgrid(coordDimY, coordDimX);

[spherCoordx, spherCoordy] = cart2pol(coordx,coordy);

mask = ones(size(coordx));

flattenCoord(:,1) = spherCoordx(mask == 1);
flattenCoord(:,2) = spherCoordy(mask == 1);

maxOrder = 2;
spherHarm = getRealSphericalHarmonics(flattenCoord, maxOrder);

spherHarm_res =reshape(spherHarm,size(mask,1),size(mask,2),size(spherHarm,2));

end

