function [ out ] = CPR_correct(data,B0,tvec,freqrange,stepsize)

B0 = B0/(2*pi); % in Hz
data(isnan(data)) = 0;

xx  = size(data,1);
yy  = size(data,2);
zz  = size(data,3);

kspace = fft32cc(data); 
        
% construct time map 
time = zeros(xx,yy);
time = repmat(tvec',[xx 1]);
time = repmat(time,[1 1 zz]);

off_res = freqrange(1):stepsize:freqrange(2);           

% use B0 map to find correct pixel in final image
count = 1;
for i = 1:xx
    for j = 1:yy
        % find which index in vec is closest to B0 per voxel i,j 
        [c k] = min(abs(off_res-B0(i,j)));
        ks(count,1) = i;
        ks(count,2) = j;
        ks(count,3) = k;
        count = count + 1;
    end
end
ks = sortrows(ks,3);

% correct the images
images = zeros(xx,yy,zz);
out    = zeros(xx,yy,zz);       
for index1 = 1:length(off_res)
    images   = ifft32cc(kspace.*exp(1i*(time)*off_res(index1)*2*pi));
    ks_index = ks(ks(:,3)==index1,:);      
    for index2 = 1:size(ks_index,1)
            out(ks_index(index2,1),ks_index(index2,2),:) = ...
                images(ks_index(index2,1),ks_index(index2,2),:);
    end
end

end

