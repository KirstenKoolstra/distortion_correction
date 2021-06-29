function [spherHarm] = getRealSphericalHarmonics(coords,maxOrder)

r0 = mean(coords(:,2));       %get the mean radius for normalisation
% spherHarm = zeros(size(coords,1),(maxOrder + 1)^2);
idx       = 1;
for n=0:maxOrder
    for m=-n:n
        if m < 0
            spherHarm(:,idx) = 1j/sqrt(2)*(coords(:,2)/r0).^n.*(harmonicY(n,m,pi/2*ones(size(coords(:,1))),coords(:,1))-((-1)^m)*harmonicY(n,-m,pi/2*ones(size(coords(:,1))),coords(:,1)));
        elseif m > 0
            spherHarm(:,idx) = 1/sqrt(2)*(coords(:,2)/r0).^n.*(harmonicY(n,-m,pi/2*ones(size(coords(:,1))),coords(:,1))+((-1)^m)*harmonicY(n,m,pi/2*ones(size(coords(:,1))),coords(:,1)));
        elseif m == 0
            spherHarm(:,idx) = (coords(:,2)/r0).^n.*harmonicY(n,m,pi/2*ones(size(coords(:,1))),coords(:,1));
        end
        idx = idx+1;
    end
end

